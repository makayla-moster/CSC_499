  #version 410											// Vertex Shader for tree.
  attribute vec3 vp;
  uniform mat4 ortho, model, view, proj;
  void main () {
    gl_Position = proj * view * model * vec4(vp, 1.0);      	// ADDING IN *ORTHO	BREAKS THE LEAVES FROM THE BRANCHES	// Tree position.
  }
