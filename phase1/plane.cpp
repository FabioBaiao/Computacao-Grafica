//SORRY FOR BAD ENGLISH :/

/*
	This is for the simplest case of the plane
	When doing the other cases, they can be used for this case and this one can be removed
*/

void plane (float l, FILE* f){

	/* 
		l must be greater than 0:
		validate here or before?
	*/

	l = l/2;

	/*
		print the following to f

	vertices:
	A -> l 0 l
	B -> l 0 -l
	C -> -l 0 -l
	D -> -l 0 l
	//number of vertices for the triangles (2 triangles each with 3 vertices)
	6
	//triangle #1
	vertex A
	vertex B
	vertex D
	//triangle #2
	vertex B
	vertex C
	vertex D
	*/
}