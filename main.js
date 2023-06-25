import {Mouse3D} from '/js/mouse3d.js';
import {GLUtil} from '/js/gl-util.js';
import {vec3, vec4, mat3, mat4} from '/js/gl-matrix-3.4.1/index.js';
import {arrayClone, getIDs, DOM, removeFromParent, hide, show, hidden} from '/js/util.js';

const ids = getIDs();
window.ids = ids;
const urlparams = new URLSearchParams(window.location.search);

let canvas;
let gl;
let glutil;
let mouse;
let wireCubeMesh;
let solidCubeMesh;
let dim = 4;
let playerStartPos = [4,4,4,-1];
let viewAngle4 = mat4.create();

//R = cos(theta) + sin(theta) *u + (1 - cos(theta)) uxu
/*
3D case:
cos(theta) -sin(theta) 0
sin(theta) cos(theta) 0
	0		0		  1

cos(theta)  0 sin(theta)
	0	    1	 0
-sin(theta) 0 cos(theta)

	1		0		0
	0	cos(theta) -sin(theta)
	0	sin(theta) cos(theta)

general case:
nz^2 + (nx^2 + ny^2) cos(theta),	-nx sin(theta),		ny sin(theta)
nx sin(theta)		,	ny^2 + (nx^2 + nz^2) cos(theta),		-nz sin(theta)
-ny sin(theta)		,		nz sin(theta)	,		nx^2 + (ny^2 + nz^2) cos(theta)

general case:
cos(theta) (xy + zx) + yz, -xy sin(theta), zx sin(theta)
xy sin(theta), (xy + yz) cos(theta) + zx, -yz sin(theta)
-zx sin(theta), yz sin(theta), (zx + yz) cos(theta) + xy

4D case:	
cos(theta) -sin(theta) 0 0
sin(theta) cos(theta)  0 0
	0         0        1 0
	0         0        0 1

cos(theta)  0 sin(theta) 0
	0       1     0      0
-sin(theta) 0 cos(theta) 0
	0       0     0      1

...

4D general case:
(xy + xz + xw) cos(theta) + yz + yw + zw, -xy sin(theta), xz sin(theta), -xw sin(theta),
xy sin(theta), (xy + yz + yw) cos(theta) + xz + xw + zw, -yz sin(theta), yw sin(theta),
-xz sin(theta), yz sin(theat), (xz + yz + zw) cos(theta) + xy + xw + yw, -zw sin(theta),
xw sin(theta), -yw sin(theta), zw sin(theta), (xw + yw + zw) cos(theta) + xy + xz + yz,

r_ij = delta_ij (a_jk^2 cos(theta) + a_kl^2 not in a_jk), sign_ij * a_ij sin(theta)
*/

function rotate4D(m, angle, xy, xz, xw, yz, yw, zw) {
	let c = Math.cos(angle);
	let s = Math.sin(angle);
	//col 1
	m[0] = c * (xy * xy + xz * xz + xw * xw) + yz * yz + yw * yw + zw * zw;
	m[1] = -s * xy;
	m[2] = s * xz;
	m[3] = -s * xw;
	//col 2
	m[4] = s * xy;
	m[5] = c * (xy * xy + yz * yz + yw * yw) + xz * xz + xw * xw + zw * zw;
	m[6] = -s * yz;
	m[7] = s * yw;
	//cal 3
	m[8] = -s * xz;
	m[9] = s * yz;
	m[10] = c * (xz * xz + yz * yz + zw * zw) + xy * xy + xw * xw + yw * yw;
	m[11] = -s * zw;
	//col 4
	m[12] = s * xw;
	m[13] = -s * yw;
	m[14] = s * zw;
	m[15] = c * (xw * xw + yw * yw + zw * zw) + xy * xy + xz * xz + yz * yz;
}

//m = matrix
//n = dimension
function orthonormalize (m, n) {
	//orthonormalize m
	for (let i = 0; i < n; ++i) {
		let len = 0;
		for (let k = 0; k < n; ++k) {
			len += m[i+n*k] * m[i+n*k];
		}
		len = 1 / Math.sqrt(len);
		for (let k = 0; k < n; ++k) {
			m[i+n*k] *= len;
		}
		//check
		for (let j = i+1; j < n; ++j) {
			let dot = 0;
			let div = 0;
			for (let k = 0; k < n; ++k) {
				dot += m[i+n*k] * m[j+n*k];
				div += m[i+n*k] * m[i+n*k];
			}
			let scale = dot / div;
			for (let k = 0; k < n; ++k) {
				m[j+n*k] -= m[i+n*k] * scale;
			}
		}
	}
}

function resize() {
	canvas.width = window.innerWidth;
	canvas.height = window.innerHeight;
	glutil.resize();

	const info = ids.info;
	const width = window.innerWidth 
		- parseInt(info.style.paddingLeft)
		- parseInt(info.style.paddingRight);
	info.style.width = width+'px';
	const height = window.innerHeight
		- parseInt(info.style.paddingTop)
		- parseInt(info.style.paddingBottom);
	info.style.height = (height - 32)+'px';
	ids.panel.style.height = (height - 16)+'px';
}

let blocks = [];
let objs = [];

class GameObject {
	constructor(args) {
		if (args.pos !== undefined) {
			this.pos = vec4.fromValues.apply(vec3, args.pos);
		} else {
			this.pos = vec4.create();
		}
	}
	draw() {
		if (!this.tex) return;
		solidCubeMesh.draw({
			uniforms : {
				pos4 : this.pos,
			},
			texs : [this.tex]
		});
	}
}

class Block extends GameObject {
}

let gravity = 9.8;
let dt = 1/20;
class DynamicObject extends GameObject {
	collisionFlags = 0;
	constructor(args) {
		super(args);
	
		if (args.vel !== undefined) {
			this.vel = vec4.fromValues.apply(vec3, args.vel);
		} else {
			this.vel = vec4.create();
		}
	}
	update() {
		if ((this.collisionFlags & (1 << 4)) == 0) {	//4 being the 1st bit of the 3rd set of bit pairs : x-,x+,y-,y+,z-,z+,w-,w+
			this.vel[2] -= gravity * dt;
		}

		vec4.scaleAndAdd(this.pos, this.pos, this.vel, dt);

		this.collisionFlags = 0;
		//check collision with all other objects ...
		// this is where chunks and bins come in handy ...
		for (let j = 0; j < blocks.length; ++j) {
			let block = blocks[j];
			let maxAbsDx = 0;
			let maxDx = undefined;
			let maxAxis = undefined;
			for (let k = dim-1; k >= 0; --k) {
				let dx = this.pos[k] - block.pos[k];
				let adx = Math.abs(dx);
				if (adx > maxAbsDx) {
					maxAbsDx = adx;
					maxDx = dx;
					maxAxis = k;
				}
			}
			if (maxAbsDx < 1) {	//colliding
				if (maxAxis === undefined) throw 'here';
				if (maxDx > 0) {
					this.pos[maxAxis] = block.pos[maxAxis] + 1;
					this.vel[maxAxis] = 0;
					this.collisionFlags |= 1 << (maxAxis << 1);	//bottom collision = 1st bit
				} else {
					this.pos[maxAxis] = block.pos[maxAxis] - 1;
					this.vel[maxAxis] = 0;
					this.collisionFlags |= 2 << (maxAxis << 1);	//top collision = 2nd bit
				}
			}
		}
	}
}

//keep the first 6 in order for bit op shortcuts
let CMD_LEFT = 1;	//x-
let CMD_RIGHT = 2;	//x+
let CMD_UP = 4;		//y-
let CMD_DOWN = 8;	//y+
let CMD_BACK = 16;	//w-
let CMD_FORTH = 32;	//w+
let CMD_MOVE_MASK = (1<<6)-1;

let CMD_JUMP = 64;	//z+

let rot3 = mat3.create();
rot3[0] = rot3[4] = rot3[9] = 1;
class Player extends DynamicObject {
	update(...args) {
		super.update(...args);

		//get 3D basis of 4D space
		mat3.fromMat4(rot3, viewAngle4);
		rot3[0+3*2] = rot3[1+3*2] = rot3[2+3*1] = rot3[2+3*0] = 0;
		rot3[2+3*2] = 1;
		mat3.transpose(rot3, rot3);
		orthonormalize(rot3, 3);

		//dirty trick: do (hyper)planar controls in xyz and vertical in w, then swap z and w
		// ... until I just outright swap the z and w in the glutil ...
		let tmp = this.vel[2];
		this.vel[2] = this.vel[3];
		this.vel[3] = tmp;
		
		this.vel[0] = this.vel[1] = this.vel[2] = 0;
		if (this.inputCmd & CMD_MOVE_MASK) {
			for (let i = 0; i < dim-1; ++i) {
				if (this.inputCmd & (1 << (i << 1))) {
					this.vel[i] = -this.speed;
				}
				if (this.inputCmd & (2 << (i << 1))) {
					this.vel[i] = this.speed;
				}
			}
			vec3.transformMat3(this.vel, this.vel, rot3);
		}
		//... and here's the swap
		tmp = this.vel[2];
		this.vel[2] = this.vel[3];
		this.vel[3] = tmp;
		
		if (this.inputCmd & CMD_JUMP) {
			if (this.collisionFlags & (1 << 4)) {
				this.vel[2] = this.jumpSpeed;
			}
		}
	
		if (this.pos[2] < -30) {
			this.reset();
		}
	}
	reset() {
		vec4.copy(this.pos, playerStartPos);
		vec4.set(this.vel, 0,0,0,0);
	}
}
Player.prototype.inputCmd = 0;
Player.prototype.speed = 2;
Player.prototype.jumpSpeed = 5;

function update() {
	//update
	for (let i = 0; i < objs.length; ++i) {
		objs[i].update();
	}
	//draw
	glutil.draw();
	for (let i = 0; i < blocks.length; ++i) {
		blocks[i].draw();
	}
	for (let i = 0; i < objs.length; ++i) {
		objs[i].draw();
	}
	window.requestAnimationFrame(update);
};

let player;
function genmap() {
	/* random crap
	for (let i = 0; i < 100; ++i) {
		let irand = (n) => { return Math.floor(n * Math.random()); };
		let pos4 = vec4.fromValues(irand(9)-4, irand(9)-4, irand(9)-4, irand(9)-4);

		blocks.push({
			pos : pos4
		});
	}
	*/
	
	/* generated level * /
	for (let i = -4; i <= 4; ++i) {
		for (let j = -4; j <= 4; ++j) {
			for (let k = -4; k <= 4; ++k) {
				for (let l = -4; l <= 4; ++l) {
					let set = true;
					set &= l <= -4;	//exclude above ground
					set &= Math.abs(i) > 1 || Math.abs(j) > 1; //exclude bottomless pit
					set |= Math.abs(i - k * Math.cos(l)) < 1 
						&& Math.abs(j - k * Math.sin(l)) < 1;
					if (set) {
						blocks.push(new Block({
							pos : [i,j,l,k]
						}));
					}
				}
			}
		}
	}
	/**/

	/* hypercube */
	let edges = [
		[-4,-4,-4],
		[4,-4,-4],
		[-4,4,-4],
		[-4,-4,4],
		[4,4,-4],
		[4,-4,4],
		[-4,4,4],
		[4,4,4]
	];
	for (let j = 0; j < 4; ++j) {
		edges.forEach(edge => {
			for (let i = -4; i <= 4; ++i) {
				let pos = arrayClone(edge);
				pos.splice(j, 0, i); 
				blocks.push(new Block({pos:pos}));
			}
		});
	}
	/**/

	/* sphere * /
	let res = 5;
	let max = 2;
	for (let i = -res; i <= res; ++i) {
		for (let j = -res; j <= res; ++j) {
			for (let k = -res; k <= res; ++k) {
				for (let l = -res; l <= res; ++l) {
					solid[i][j][k][l] = i*i + j*j + k*k + l*l < res*res - 1;
				}
			}
		}
	}

	/* mandel-julia fractal * /
	let res = 5;
	let max = 2;
	let maxiter = 10;
	let solid = [];
	for (let i = -res; i <= res; ++i) {
		solid[i] = [];
		for (let j = -res; j <= res; ++j) {
			solid[i][j] = [];
			for (let k = -res; k <= res; ++k) {
				solid[i][j][k] = [];
				for (let l = -res; l <= res; ++l) {
					let zr = i*max/res;
					let zi = j*max/res;
					let zj = k*max/res;
					let zk = l*max/res;
					let cr = -1;
					let ci = .2;
					let cj = 0;
					let ck = 0;
					let iter = 0;
					for (; iter < maxiter; ++iter) {
						let zrn = zr * zr - zi * zi + cr;
						let zin = 2 * zr * zi + ci;
						zr = zrn;
						zi = zin;
						if (zr * zr + zi * zi > 4) break;
					}
					solid[i][j][k][l] = iter == maxiter;
				}
			}
		}
	}
	let issolid = (i,j,k,l) => {
		let r = solid;
		if (!(i in r)) return; r = r[i];
		if (!(j in r)) return; r = r[j];
		if (!(k in r)) return; r = r[k];
		if (!(l in r)) return; r = r[l];
		return r;
	};
	for (let i = -res; i <= res; ++i) {
		for (let j = -res; j <= res; ++j) {
			for (let k = -res; k <= res; ++k) {
				for (let l = -res; l <= res; ++l) {
					if (issolid(i,j,k,l)) {
						if (issolid(i+1,j,k,l) &&
							issolid(i-1,j,k,l) &&
							issolid(i,j+1,k,l) &&
							issolid(i,j-1,k,l) &&
							issolid(i,j,k+1,l) &&
							issolid(i,j,k-1,l) &&
							issolid(i,j,k,l+1) &&
							issolid(i,j,k,l-1))
						{
							blocks.push(new Block({pos:[i,j,k,l]}));
						}
					}
				}
			}
		}
	}
	/**/
/*
	objs.push(player = new Player({
		pos : playerStartPos
	}));
*/
}

function onkeydown(e) {
	if (!player) return;
	let keyCode = e.keyCode;
	switch (keyCode) {
	case 87:	//'w'
	case 38:	//up
		player.inputCmd |= CMD_DOWN;
		break;
	case 83:	//'s'
	case 40:	//down
		player.inputCmd |= CMD_UP;
		break;
	case 65:	//'a'
	case 37:	//left
		player.inputCmd |= CMD_LEFT;
		break;
	case 68:	//'d'
	case 39:	//right
		player.inputCmd |= CMD_RIGHT;
		break;
	case 81:	//'q'
		player.inputCmd |= CMD_BACK;
		break;
	case 69:	//'e'
		player.inputCmd |= CMD_FORTH;
		break;
	case 13:	//enter
	case 32:	//space
		player.inputCmd |= CMD_JUMP;
		break;
	default:
		return;	//...and don't prevent default
	}
	e.preventDefault();
}

function onkeyup(e) {
	if (!player) return;
	let keyCode = e.keyCode;
	switch (keyCode) {
	case 87:	//'w'
	case 38:	//up
		player.inputCmd &= ~CMD_DOWN;
		break;
	case 83:	//'s'
	case 40:	//down
		player.inputCmd &= ~CMD_UP;
		break;
	case 65:	//'a'
	case 37:	//left
		player.inputCmd &= ~CMD_LEFT;
		break;
	case 68:	//'d'
	case 39:	//right
		player.inputCmd &= ~CMD_RIGHT;
		break;
	case 13:	//enter
	case 32:	//space
		player.inputCmd &= ~CMD_JUMP;
		break;
	case 81:	//'q'
		player.inputCmd &= ~CMD_BACK;
		break;
	case 69:	//'e'
		player.inputCmd &= ~CMD_FORTH;
		break;
	default:
		return;	//...and don't prevent default
	}
	e.preventDefault();
}

function onkeypress(e) {
	e.preventDefault();
}

function initInput() {
	window.addEventListener('keydown', onkeydown);
	window.addEventListener('keyup', onkeyup);
	window.addEventListener('keypress', onkeypress);
}



ids.panelButton.addEventListener('click', e => {
	if (hidden(ids.panel)) {
		show(ids.panel);
		hide(ids.info);
	} else {
		hide(ids.panel);
	}
});
ids.infoButton.addEventListener('click', e => {
	if (hidden(ids.info)) {
		show(ids.info);
		hide(ids.panel);
	} else {
		hide(ids.info);
	}
});

canvas = DOM(
	'canvas', {
	css : {
		left : 0,
		top : 0,
		position : 'absolute',
		userSelect : 'none',
	},
	prependTo : document.body,
});
window.canvas = canvas;
try {
	glutil = new GLUtil({canvas:canvas});
	gl = glutil.context;
} catch (e) {
	removeFromParent(canvas);
	show(ids.webglfail);
	throw e;
}
show(ids.menu);
show(ids.panel);

if (urlparams.get('info')) {
	show(ids.info);
	hide(ids.panel);
}

/*
ids.resetPlayer.addEventListener('click', e => {
	if (!player) return;
	player.reset();
});
*/

let rotationMethod = 'xyw';
const rotation = document.getElementById('rotation');
rotation.addEventListener('change', e => {
	rotationMethod = rotation.value;
});

let originalViewAngle4 = [
	1, 0, 0, 0, 
	0, 0.5, -0.866, 0, 
	0, 0.866, 0.5, 0, 
	0, 0, 0, 1
];
ids.resetRotation.addEventListener('click', e => {
	mat4.copy(viewAngle4, originalViewAngle4);
});
mat4.copy(viewAngle4, originalViewAngle4);
orthonormalize(viewAngle4, 4);

let viewPos4 = vec4.create();

glutil.view.zNear = .1;
glutil.view.zFar = 100;
glutil.onfps = (fps) => { ids.fps.innerText = fps; };
viewPos4[2] = 10;

gl.enable(gl.DEPTH_TEST);
gl.clearColor(.6, .8, 1., 1.);

const tesseractVertexCode = 
`in vec4 vertex;
in vec2 texCoord;
uniform mat4 projMat;
uniform vec4 pos4;
uniform vec4 viewPos4;
uniform mat4 viewAngle4;
out vec4 srcVertex;
out vec2 texCoordV;
void main() {
	texCoordV = texCoord;
	srcVertex = vertex;
	vec4 vtx4 = vertex;
	vtx4 += pos4;
	vtx4 = viewAngle4 * vtx4;
	vtx4 -= viewPos4;
	vtx4.xyz *= 1. / (1. + vtx4.w * vtx4.w);
	vtx4.w = 1.;
	gl_Position = projMat * vtx4;
}
`;

//application of w locally.  good if were not rotating in the w plane
const slicesVertexCode = 
`in vec4 vertex;
in vec2 texCoord;
uniform mat4 projMat;
uniform vec4 pos4;
uniform vec4 viewPos4;
uniform mat4 viewAngle4;
out vec4 srcVertex;
out vec2 texCoordV;

highp mat4 mat4orthonormal(highp mat4 m) {
	m[0] /= length(m[0]);
	m[1] -= (dot(m[0], m[1]) / dot(m[0], m[0])) * m[0];
	m[2] -= (dot(m[0], m[2]) / dot(m[0], m[0])) * m[0];
	m[3] -= (dot(m[0], m[3]) / dot(m[0], m[0])) * m[0];
	m[1] /= length(m[1]);
	m[2] -= (dot(m[1], m[2]) / dot(m[1], m[1])) * m[1];
	m[3] -= (dot(m[1], m[3]) / dot(m[1], m[1])) * m[1];
	m[2] /= length(m[2]);
	m[3] -= (dot(m[2], m[3]) / dot(m[2], m[2])) * m[2];
	m[3] /= length(m[3]);
	return m;
}

void main() {
	texCoordV = texCoord;
	srcVertex = vertex;
	vec4 vtx4 = vertex;

	mat4 rot3 = viewAngle4;
	rot3[0][3] = 0.;
	rot3[1][3] = 0.;
	rot3[2][3] = 0.;
	rot3[3][0] = 0.;
	rot3[3][1] = 0.;
	rot3[3][2] = 0.;
	rot3[3][3] = 1.;
	rot3 = mat4orthonormal(rot3);

	mat4 rot4 = viewAngle4;
	rot4[0].xyz = vec3(1., 0., 0.);
	rot4[1].xyz = vec3(0., 1., 0.);
	rot4[2].xyz = vec3(0., 0., 1.);
	rot4 = mat4orthonormal(rot4);

	//apply w-component of object before applying w-component of vertex
	//apply w-component of rotation here too
	vec4 pos4r = rot4 * pos4;//viewAngle4 * pos4;

	//grid
	float f = 3. * floor((pos4r.w + 1.) / 3.);	//-3, 0, 3, closest div 3 rounded down ...
	float g = 3. * floor((pos4r.w + 1.) - f);	//0, 1, 2, 3
	vtx4.x += 3.5 * f;
	vtx4.y += 3.5 * (g - 3.);

	vtx4.xyz += pos4r.xyz;

	//apply 3D component of rotation
	vtx4 = rot3 * vtx4;

	vtx4 -= viewPos4;
	vtx4.xyz *= 1. / (1. + vtx4.w * vtx4.w);
	vtx4.w = 1.;
	gl_Position = projMat * vtx4;
}
`;

const fragmentCode = 
`uniform sampler2D tex;
in vec4 srcVertex;
in vec2 texCoordV;
out vec4 fragColor;
void main() {
	fragColor = vec4(0.);
	//hue by dimension
	fragColor.r += .5 + srcVertex.w;
	fragColor.b += .5 - srcVertex.w;
	fragColor.r += .5 + srcVertex.x;
	fragColor.g += .5 - srcVertex.x;
	fragColor.b += .5 + srcVertex.y;
	fragColor.r += .5 - srcVertex.y;
	fragColor.g += .5 + srcVertex.z;
	fragColor.b += .5 - srcVertex.z;
	fragColor *= .1;	
	
	//apply base color
	fragColor += texture(tex, texCoordV);
	
	//for kicks
	fragColor.a = 1.;

	float light = dot(normalize(srcVertex.xyz), vec3(1.));
	light = max(light, .5);
	fragColor *= light;
}
`;

let tesseractShader = new glutil.Program({
	vertexCode : tesseractVertexCode,	
	fragmentCode : fragmentCode,
});
	
let slicesShader = new glutil.Program({
	vertexCode : slicesVertexCode,
	fragmentCode : fragmentCode,
});

let scale = 1;//.95;
/*
let vertex = [];
for (let i = 0; i < 1<<dim; ++i) {
	for (let j = 0; j < dim; ++j) {
		vertex.push((((i>>j)&1)-.5) * scale);
	}
}
let edges = [];	//indexes for edges
for (let i = 0; i < 1<<dim; ++i) {
	for (let j = 0; j < dim; ++j) {
		let ea = i;
		let eb = i ^ (1 << j);
		if (eb > ea) {
			edges.push(ea);
			edges.push(eb);
		}
	}
}

wireCubeMesh = new glutil.SceneObject({
	mode : gl.LINES,
	shader : shader,
	indexes : new glutil.ElementArrayBuffer({
		data : edges
	}),
	attrs : {
		vertex : new glutil.ArrayBuffer({
			dim : dim,
			data : vertex
		})
	},
	uniforms : {
		viewPos4 : viewPos4,		//view position
		viewAngle4 : viewAngle4		//view rotation
	},
	static : true 
});
*/

let quadVertexes = [];	//indexes for quads ... er, triangles
let texCoords = [];
for (let i = 0; i < 1<<dim; ++i) {
	for (let j = 0; j < dim-1; ++j) {
		for (let k = j+1; k < dim; ++k) {
			let q1 = i;
			let q2 = i ^ (1 << j);
			let q3 = i ^ (1 << j) ^ (1 << k);
			let q4 = i ^ (1 << k);
			if (q3 > q2 && q2 > q1) {
				let push = (n) => {
					for (let k = 0; k < dim; ++k) {
						quadVertexes.push((((n>>k)&1)-.5) * scale);
					}
				}
				push(q1);
				push(q2);
				push(q3);
				push(q3);
				push(q4);
				push(q1);
				texCoords = texCoords.concat([0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0]);
			}
		}
	}
}

new glutil.Texture2D({
	url : 'tex/bricks.png',
	minFilter : gl.NEAREST,
	magFilter : gl.LINEAR,
	onload : function() { Block.prototype.tex = this; },
});

new glutil.Texture2D({
	url : 'tex/player.png',
	minFilter : gl.NEAREST,
	magFilter : gl.LINEAR,
	onload : function() { Player.prototype.tex = this; },
});

ids.renderer.addEventListener('change', e => {
	solidCubeMesh.shader = ({
		tesseract : tesseractShader,
		slices : slicesShader
	})[ids.renderer.value];
});

solidCubeMesh = new glutil.SceneObject({
	mode : gl.TRIANGLES,
	shader : tesseractShader,
	attrs : {
		vertex : new glutil.ArrayBuffer({
			dim : dim,
			data : quadVertexes,
		}),
		texCoord : new glutil.ArrayBuffer({
			dim : 2, 
			data : texCoords
		})
	},
	uniforms : {
		viewPos4 : viewPos4,
		viewAngle4 : viewAngle4,
		tex : 0,
	},
	texs : [],
	static : true
});

genmap();

let tmpR = mat4.create();
mouse = new Mouse3D({
	pressObj : canvas,
	move : (dx,dy) => {
		let rotAngle = Math.PI / 180 * Math.sqrt(dx*dx + dy*dy);
		let r = Math.sqrt(dx*dx + dy*dy);
		if (r == 0) return;
		if (rotationMethod == '4d') {
			rotate4D(tmpR, -.01 * dx, 0, 0, 1, 0, 0, 0);	//xw rotation
			mat4.mul(viewAngle4, viewAngle4, tmpR);
			rotate4D(tmpR, -.01 * dy, 0, 0, 0, 0, 1, 0);	//yw rotation
			mat4.mul(viewAngle4, viewAngle4, tmpR);
			orthonormalize(viewAngle4, 4);
		} else if (rotationMethod == '3d') {
			rotate4D(tmpR, -.01 * dx, 1, 0, 0, 0, 0, 0);	//xy rotation
			mat4.mul(viewAngle4, viewAngle4, tmpR);
			rotate4D(tmpR, -.01 * dy, 0, 0, 0, 1, 0, 0);	//yz rotation
			mat4.mul(viewAngle4, tmpR, viewAngle4);
			orthonormalize(viewAngle4, 4);
		} else {
			let xy = 0;
			let xz = 0;
			let xw = 0;
			let yz = 0;
			let yw = 0;
			let zw = 0;

			if (rotationMethod == 'xyz') {
				xz = -dx / r;
				yz = -dy / r;
			} else if (rotationMethod == 'xyw') {
				xw = -dx / r;
				yw = -dy / r;
			} else if (rotationMethod == 'xzw') {
				xw = -dx / r;
				zw = -dy / r;
			} else if (rotationMethod == 'yzw') {
				zw = -dx / r;
				yw = -dy / r;
			}

			rotate4D(tmpR, rotAngle, xy, xz, xw, yz, yw, zw);
			mat4.mul(viewAngle4, tmpR, viewAngle4);
			orthonormalize(viewAngle4, 4);
		}
	},
	zoom : (dz) => {
		viewPos4[2] += .001 * dz;
	}
});

initInput();

window.addEventListener('resize', resize);
resize();

update();
