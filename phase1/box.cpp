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

void box2 (float x, float y, float z, float nDiv){

	/*
		x, y, z and nDiv must be greater than 0.
	*/

	float deltaX = x / nDiv;
	float deltaY = y / nDiv;
	float deltaZ = z / nDiv;

	float xi, yi, zi;

	x = x/2;
	y = y/2;
	z = z/2;

	/*
		FRONT FACE

		The xi, yi, zi are the initial value to print each rectangle (not each face).
		Their starting value is the bottom right corner of the face to be printed.
	*/
	xi = x;
	yi = -y;
	zi = z;
	/* 
		The outer cycle iterates through each line of the face (down -> up).
	*/
	for (int i = 0; i < nDiv; i++){
		/*
			The inner cycle prints all rectangles of the same line in the face being print (from the rigth to the left).
		*/
		for (int j = 0; j < nDiv; j++){
			/*
				xi yi zi; //(starting vertex)
				xi (yi+deltaY) zi; //(above vertex)
				(xi-deltaX) (yi+deltaY) zi; //(left-up diagonal vertex)

				xi yi zi; //(starting vertex)
				(xi-deltaX) (yi+deltaY) zi; //(left-up diagonal vertex)
				(xi-deltaX) yi zi; //(left vertex)
			*/
			xi -= deltaX; // the next iteration is at the left of the starting value.
		}
		/*
			To start a new line, it's necessary to go back to the inital bottom right corner but one line up. 
		*/
		xi += nDiv * deltaX;
		yi += deltaY;
	}
	/*
		To make the other faces the process is similar, but the coordinates to sum and subtract are different.

		IF YOU AGREE WITH THIS ALGORITHM, I (Fabio Baiao) VOLUNTEER (as tribute :D) TO COMPLETE THE PSEUDOCODE OF THE OTHER FACES.
	*/
}