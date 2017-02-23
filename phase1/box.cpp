//AGAIN, SORRY FOR BAD ENGLISH

/*
	This pseudocode is for the box without division numbers.
	When adding the code for the case with division numbers, we can remove this one, and use the new with division number = 1
*/

void box (float x, float y, float z, FILE *f){
	
	/*
		x, y and z must be greater than 0
	*/

	x = x/2;
	y = y/2;
	z = z/2;

	/*

	I'm assuming that the given coordinates are according to the OpenGL coordinate system
	vertices:
	A -> x -y z
	B -> x y z
	C -> -x y z
	D -> -x -y z
	E -> x -y -z
	F -> x y -z
	G -> -x y -z
	H -> -x -y -z
	//number of vertices for the triangles (6 faces each with 2 triangles each with 3 vertices)
	36
	//front face
	//triangle #1
	vertex A
	vertex B
	vertex D
	//triangle #2
	vertex B
	vertex C
	vertex D
	//back face
	//triangle #1
	vertex E
	vertex H
	vertex F
	//triangle #2
	vertex F
	vertex H
	vertex G
	//up face
	//trangle #1
	vertex B
	vertex F
	vertex C
	//triangle #2
	vertex C
	vertex F
	vertex G
	//down face
	//triangle #1
	vertex A
	vertex D
	vertex E
	//triangle #2
	vertex D
	vertex H
	vertex E
	//left face
	//triangle #1
	vertex C
	vertex H
	vertex D
	//triangle #2
	vertex C
	vertex G
	vertex H
	//right face
	//triangle #1
	vertex A
	vertex E
	vertex B
	//triangle #2
	vertex B
	vertex E
	vertex F
	*/
}