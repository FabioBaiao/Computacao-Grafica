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

	//number of vertices
	6
	//triangle #1
	l 0 l (vertex A)
	l 0 -l (vertex B)
	-l 0 l (vertex D)
	//triangle #2
	l 0 -l (vertex B)
	-l 0 -l (vertex C)
	-l 0 l (vertex D)
	*/
}