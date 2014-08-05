var canvas;
var gl;
var renderer;
var mouse;
var wireCubeMesh;
var solidCubeMesh;
var dim = 4;
var playerStartPos = [4,4,4,-1];
var viewAngle4 = mat4.create();

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

mat4.rotate4D = function(m, angle, xy, xz, xw, yz, yw, zw) {
	var c = Math.cos(angle);
	var s = Math.sin(angle);
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
orthonormalize = function(m, n) {
	//orthonormalize m
	for (var i = 0; i < n; ++i) {
		var len = 0;
		for (var k = 0; k < n; ++k) {
			len += m[i+n*k] * m[i+n*k];
		}
		len = 1 / Math.sqrt(len);
		for (var k = 0; k < n; ++k) {
			m[i+n*k] *= len;
		}
		//check
		for (var j = i+1; j < n; ++j) {
			var dot = 0;
			var div = 0;
			for (var k = 0; k < n; ++k) {
				dot += m[i+n*k] * m[j+n*k];
				div += m[i+n*k] * m[i+n*k];
			}
			var scale = dot / div;
			for (var k = 0; k < n; ++k) {
				m[j+n*k] -= m[i+n*k] * scale;
			}
		}
	}
}

function resize() {
	canvas.width = window.innerWidth;
	canvas.height = window.innerHeight;
	renderer.resize();

	var info = $('#info');
	var width = window.innerWidth 
		- parseInt(info.css('padding-left'))
		- parseInt(info.css('padding-right'));
	info.width(width);
	var height = window.innerHeight
		- parseInt(info.css('padding-top'))
		- parseInt(info.css('padding-bottom'));
	info.height(height - 32);
	$('#panel').height(height - 16);
}

var blocks = [];
var objs = [];

var GameObject = makeClass({
	init : function(args) {
		if (args.pos !== undefined) {
			this.pos = vec4.fromValues.apply(vec3, args.pos);
		} else {
			this.pos = vec4.create();
		}
	},
	draw : function() {
		solidCubeMesh.draw({
			uniforms : {
				pos4 : this.pos
			},
			texs : [this.tex]
		});
	}
});

var Block = makeClass({
	super : GameObject
});

var gravity = 9.8;
var dt = 1/20;
var DynamicObject = makeClass({
	super : GameObject,
	collisionFlags : 0,
	init : function(args) {
		DynamicObject.super.apply(this, arguments);
	
		if (args.vel !== undefined) {
			this.vel = vec4.fromValues.apply(vec3, args.vel);
		} else {
			this.vel = vec4.create();
		}
	},
	update : function() {
		if ((this.collisionFlags & (1 << 4)) == 0) {	//4 being the 1st bit of the 3rd set of bit pairs : x-,x+,y-,y+,z-,z+,w-,w+
			this.vel[2] -= gravity * dt;
		}

		vec4.scaleAndAdd(this.pos, this.pos, this.vel, dt);

		this.collisionFlags = 0;
		//check collision with all other objects ...
		// this is where chunks and bins come in handy ...
		for (var j = 0; j < blocks.length; ++j) {
			var block = blocks[j];
			var maxAbsDx = 0;
			var maxDx = undefined;
			var maxAxis = undefined;
			for (var k = dim-1; k >= 0; --k) {
				var dx = this.pos[k] - block.pos[k];
				var adx = Math.abs(dx);
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
});

//keep the first 6 in order for bit op shortcuts
var CMD_LEFT = 1;	//x-
var CMD_RIGHT = 2;	//x+
var CMD_UP = 4;		//y-
var CMD_DOWN = 8;	//y+
var CMD_BACK = 16;	//w-
var CMD_FORTH = 32;	//w+
var CMD_MOVE_MASK = (1<<6)-1;

var CMD_JUMP = 64;	//z+

var rot3 = mat3.create();
rot3[0] = rot3[4] = rot3[9] = 1;
var Player = makeClass({
	super : DynamicObject,
	inputCmd : 0,
	speed : 2,
	jumpSpeed : 5,
	update : function() {
		Player.superProto.update.apply(this, arguments);

		//get 3D basis of 4D space
		mat3.fromMat4(rot3, viewAngle4);
		rot3[0+3*2] = rot3[1+3*2] = rot3[2+3*1] = rot3[2+3*0] = 0;
		rot3[2+3*2] = 1;
		mat3.transpose(rot3, rot3);
		orthonormalize(rot3, 3);

		//dirty trick: do (hyper)planar controls in xyz and vertical in w, then swap z and w
		// ... until I just outright swap the z and w in the renderer ...
		var tmp = this.vel[2];
		this.vel[2] = this.vel[3];
		this.vel[3] = tmp;
		
		this.vel[0] = this.vel[1] = this.vel[2] = 0;
		if (this.inputCmd & CMD_MOVE_MASK) {
			for (var i = 0; i < dim-1; ++i) {
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
		var tmp = this.vel[2];
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
	},
	reset : function() {
		vec4.copy(this.pos, playerStartPos);
		vec4.set(this.vel, 0,0,0,0);
	}
});


function update() {
	//update
	for (var i = 0; i < objs.length; ++i) {
		objs[i].update();
	}
	//draw
	renderer.draw();
	for (var i = 0; i < blocks.length; ++i) {
		blocks[i].draw();
	}
	for (var i = 0; i < objs.length; ++i) {
		objs[i].draw();
	}
	requestAnimFrame(update);
};

var player;
function genmap() {
	/* random crap
	for (var i = 0; i < 100; ++i) {
		var irand = function(n) { return Math.floor(n * Math.random()); };
		var pos4 = vec4.fromValues(irand(9)-4, irand(9)-4, irand(9)-4, irand(9)-4);

		blocks.push({
			pos : pos4
		});
	}
	*/
	
	/* generated level * /
	for (var i = -4; i <= 4; ++i) {
		for (var j = -4; j <= 4; ++j) {
			for (var k = -4; k <= 4; ++k) {
				for (var l = -4; l <= 4; ++l) {
					var set = true;
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
	var edges = [
		[-4,-4,-4],
		[4,-4,-4],
		[-4,4,-4],
		[-4,-4,4],
		[4,4,-4],
		[4,-4,4],
		[-4,4,4],
		[4,4,4]
	];
	for (var j = 0; j < 4; ++j) {
		$.each(edges, function(_,edge) {
			for (var i = -4; i <= 4; ++i) {
				var pos = Array.apply(undefined, edge);
				pos.splice(j, 0, i); 
				blocks.push(new Block({pos:pos}));
			}
		});
	}
	/**/

	/* sphere * /
	var res = 5;
	var max = 2;
	for (var i = -res; i <= res; ++i) {
		for (var j = -res; j <= res; ++j) {
			for (var k = -res; k <= res; ++k) {
				for (var l = -res; l <= res; ++l) {
					solid[i][j][k][l] = i*i + j*j + k*k + l*l < res*res - 1;
				}
			}
		}
	}

	/* mandel-julia fractal * /
	var res = 5;
	var max = 2;
	var maxiter = 10;
	var solid = [];
	for (var i = -res; i <= res; ++i) {
		solid[i] = [];
		for (var j = -res; j <= res; ++j) {
			solid[i][j] = [];
			for (var k = -res; k <= res; ++k) {
				solid[i][j][k] = [];
				for (var l = -res; l <= res; ++l) {
					var zr = i*max/res;
					var zi = j*max/res;
					var zj = k*max/res;
					var zk = l*max/res;
					var cr = -1;
					var ci = .2;
					var cj = 0;
					var ck = 0;
					var iter = 0;
					for (; iter < maxiter; ++iter) {
						var zrn = zr * zr - zi * zi + cr;
						var zin = 2 * zr * zi + ci;
						zr = zrn;
						zi = zin;
						if (zr * zr + zi * zi > 4) break;
					}
					solid[i][j][k][l] = iter == maxiter;
				}
			}
		}
	}
	var issolid = function(i,j,k,l) {
		var r = solid;
		if (!(i in r)) return; r = r[i];
		if (!(j in r)) return; r = r[j];
		if (!(k in r)) return; r = r[k];
		if (!(l in r)) return; r = r[l];
		return r;
	};
	for (var i = -res; i <= res; ++i) {
		for (var j = -res; j <= res; ++j) {
			for (var k = -res; k <= res; ++k) {
				for (var l = -res; l <= res; ++l) {
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

function onkeydown(event) {
	if (!player) return;
	var keyCode = event.keyCode;
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
	event.preventDefault();
}

function onkeyup(event) {
	if (!player) return;
	var keyCode = event.keyCode;
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
	event.preventDefault();
}

function onkeypress(event) {
	event.preventDefault();
}

function initInput() {
	$(window)
		.keydown(onkeydown)
		.keyup(onkeyup)
		.keypress(onkeypress);
}

$(document).ready(function() {
	$('#panelButton').click(function() {
		var panel = $('#panel');	
		if (panel.css('display') == 'none') {
			panel.show();
			$('#info').hide();
		} else {
			panel.hide();
		}
	});
	$('#infoButton').click(function() {
		var info = $('#info');
		if (info.css('display') == 'none') {
			info.show();
			$('#panel').hide();
		} else {
			info.hide();
		}
	});
	
	canvas = $('<canvas>', {
		css : {
			left : 0,
			top : 0,
			position : 'absolute'
		}
	}).prependTo(document.body).get(0);
	$(window).disableSelection()

	try {
		renderer = new GL.CanvasRenderer({canvas:canvas});
		gl = renderer.gl;
	} catch (e) {
		$(canvas).remove();
		$('#webglfail').show();
		throw e;
	}
	$('#menu').show();
	$('#panel').show();

	if ($.url().param('info')) {
		$('#info').show();
		$('#panel').hide();
	}

	/*
	$('#resetPlayer').click(function() {
		if (!player) return;
		player.reset();
	});
	*/
	
	rotationMethod = '3d';
	$('#rotation').change(function() {
		rotationMethod = $(this).val();
	});
	
	var originalViewAngle4 = [
		1, 0, 0, 0, 
		0, 0.5, -0.866, 0, 
		0, 0.866, 0.5, 0, 
		0, 0, 0, 1
	];
	$('#resetRotation').click(function() {
		mat4.copy(viewAngle4, originalViewAngle4);
	});
	mat4.copy(viewAngle4, originalViewAngle4);
	orthonormalize(viewAngle4, 4);
	
	viewPos4 = vec4.create();
	
	renderer.view.zNear = .1;
	renderer.view.zFar = 100;
	renderer.onfps = function(fps) { $('#fps').text(fps); };
	viewPos4[2] = 10;

	gl.enable(gl.DEPTH_TEST);
	gl.clearColor(.6, .8, 1., 1.);


	var tesseractVertexShader = new GL.VertexShader({
		code : GL.vertexPrecision + mlstr(function(){/*
attribute vec4 vertex;
attribute vec2 texCoord;
uniform mat4 projMat;
uniform vec4 pos4;
uniform vec4 viewPos4;
uniform mat4 viewAngle4;
varying vec4 srcVertex;
varying vec2 texCoordV;
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
*/})
	});

	//application of w locally.  good if were not rotating in the w plane
	var slicesVertexShader = new GL.VertexShader({
		code : GL.vertexPrecision + mlstr(function(){/*
attribute vec4 vertex;
attribute vec2 texCoord;
uniform mat4 projMat;
uniform vec4 pos4;
uniform vec4 viewPos4;
uniform mat4 viewAngle4;
varying vec4 srcVertex;
varying vec2 texCoordV;

//http://stackoverflow.com/questions/18034677/transpose-a-mat4-in-opengl-es-2-0-glsl
highp mat4 transpose(in highp mat4 inMatrix) {
	highp vec4 i0 = inMatrix[0];
	highp vec4 i1 = inMatrix[1];
	highp vec4 i2 = inMatrix[2];
	highp vec4 i3 = inMatrix[3];

	highp mat4 outMatrix = mat4(
                 vec4(i0.x, i1.x, i2.x, i3.x),
                 vec4(i0.y, i1.y, i3.y, i3.y),
                 vec4(i0.z, i1.z, i3.z, i3.z),
                 vec4(i0.w, i1.w, i3.w, i3.w)
                 );
	return outMatrix;
}

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
*/})
	});

	var fragmentShader = new GL.FragmentShader({
		code : GL.fragmentPrecision + mlstr(function(){/*
varying vec4 srcVertex;
varying vec2 texCoordV;
uniform sampler2D tex;
void main() {
	gl_FragColor = vec4(0.);
	//hue by dimension
	gl_FragColor.r += .5 + srcVertex.w;
	gl_FragColor.b += .5 - srcVertex.w;
	gl_FragColor.r += .5 + srcVertex.x;
	gl_FragColor.g += .5 - srcVertex.x;
	gl_FragColor.b += .5 + srcVertex.y;
	gl_FragColor.r += .5 - srcVertex.y;
	gl_FragColor.g += .5 + srcVertex.z;
	gl_FragColor.b += .5 - srcVertex.z;
	gl_FragColor *= .1;	
	
	//apply base color
	gl_FragColor += texture2D(tex, texCoordV);
	
	//for kicks
	gl_FragColor.a = 1.;

	float light = dot(normalize(srcVertex.xyz), vec3(1.));
	light = max(light, .5);
	gl_FragColor *= light;
}
*/})
	});

	var tesseractShader = new GL.ShaderProgram({
		vertexShader : tesseractVertexShader,	
		fragmentShader : fragmentShader
	});
		
	var slicesShader = new GL.ShaderProgram({
		vertexShader : slicesVertexShader,
		fragmentShader : fragmentShader
	});

	var scale = 1;//.95;
/*
	var vertex = [];
	for (var i = 0; i < 1<<dim; ++i) {
		for (var j = 0; j < dim; ++j) {
			vertex.push((((i>>j)&1)-.5) * scale);
		}
	}
	var edges = [];	//indexes for edges
	for (var i = 0; i < 1<<dim; ++i) {
		for (var j = 0; j < dim; ++j) {
			var ea = i;
			var eb = i ^ (1 << j);
			if (eb > ea) {
				edges.push(ea);
				edges.push(eb);
			}
		}
	}

	wireCubeMesh = new GL.SceneObject({
		mode : gl.LINES,
		shader : shader,
		indexes : new GL.ElementArrayBuffer({
			data : edges
		}),
		attrs : {
			vertex : new GL.ArrayBuffer({
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

	var quadVertexes = [];	//indexes for quads ... er, triangles
	var texCoords = [];
	for (var i = 0; i < 1<<dim; ++i) {
		for (var j = 0; j < dim-1; ++j) {
			for (var k = j+1; k < dim; ++k) {
				var q1 = i;
				var q2 = i ^ (1 << j);
				var q3 = i ^ (1 << j) ^ (1 << k);
				var q4 = i ^ (1 << k);
				if (q3 > q2 && q2 > q1) {
					var push = function(n) {
						for (var k = 0; k < dim; ++k) {
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
	
	new GL.Texture2D({
		url : 'tex/bricks.png',
		minFilter : gl.NEAREST,
		magFilter : gl.LINEAR,
		onload : function() {
			Block.prototype.tex = this;
		}
	});

	new GL.Texture2D({
		url : 'tex/player.png',
		minFilter : gl.NEAREST,
		magFilter : gl.LINEAR,
		onload : function() {
			Player.prototype.tex = this;
		}
	});

	$('#renderer').change(function() {
		solidCubeMesh.shader = ({
			tesseract : tesseractShader,
			slices : slicesShader
		})[$(this).val()];
	});

	solidCubeMesh = new GL.SceneObject({
		mode : gl.TRIANGLES,
		shader : tesseractShader,
		attrs : {
			vertex : new GL.ArrayBuffer({dim : dim, data : quadVertexes, keep : true}),
			texCoord : new GL.ArrayBuffer({dim : 2, data : texCoords, keep : true})
		},
		uniforms : {
			viewPos4 : viewPos4,
			viewAngle4 : viewAngle4,
			tex : 0
		},
		texs : [],
		static : true
	});

	genmap();

	var tmpR = mat4.create();
	mouse = new Mouse3D({
		pressObj : canvas,
		move : function(dx,dy) {
			var rotAngle = Math.PI / 180 * Math.sqrt(dx*dx + dy*dy);
			var r = Math.sqrt(dx*dx + dy*dy);
			if (r == 0) return;
			if (rotationMethod == '4d') {
				mat4.rotate4D(tmpR, -.01 * dx, 0, 0, 1, 0, 0, 0);	//xw rotation
				mat4.mul(viewAngle4, viewAngle4, tmpR);
				mat4.rotate4D(tmpR, -.01 * dy, 0, 0, 0, 0, 1, 0);	//yw rotation
				mat4.mul(viewAngle4, viewAngle4, tmpR);
				orthonormalize(viewAngle4, 4);
			} else if (rotationMethod == '3d') {
				mat4.rotate4D(tmpR, -.01 * dx, 1, 0, 0, 0, 0, 0);	//xy rotation
				mat4.mul(viewAngle4, viewAngle4, tmpR);
				mat4.rotate4D(tmpR, -.01 * dy, 0, 0, 0, 1, 0, 0);	//yz rotation
				mat4.mul(viewAngle4, tmpR, viewAngle4);
				orthonormalize(viewAngle4, 4);
			} else {
				var xy = 0;
				var xz = 0;
				var xw = 0;
				var yz = 0;
				var yw = 0;
				var zw = 0;

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

				mat4.rotate4D(tmpR, rotAngle, xy, xz, xw, yz, yw, zw);
				mat4.mul(viewAngle4, tmpR, viewAngle4);
				orthonormalize(viewAngle4, 4);
			}
		},
		zoom : function(dz) {
			viewPos4[2] += .001 * dz;
		}
	});

	initInput();

	$(window).resize(resize);
	resize();
	
	update();
});

