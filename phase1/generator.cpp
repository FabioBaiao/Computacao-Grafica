#include <iostream>
#include <fstream>
#include <regex>

using namespace std;

bool validate_input(int, char**);
void generate_sphere(char**);
void generate_plane(float, float, float,char*);

std::string drawSphere(float radius, float slices, int stacks) {

	std:string ret = "";
	int count = 0;
	int half_stacks = (int) stacks / 2;
	float b_step = M_PI / stacks;
	float a_step = 2*M_PI / slices;
		
	for(float i = M_PI/2; i > -M_PI/2; i -= b_step)
	{
//		float beta = i * (M_PI) / stacks;
		float delta_b = i - (M_PI) / stacks;
		
		for(int j = 0; j < slices; j++)
		{
			float alfa = j * (2*M_PI)/slices;
			float delta_a = alfa + (2*M_PI)/slices;

			ret += std::to_string(radius*cos(i)*sin(alfa)) + " " + 
				   std::to_string(radius*sin(i)) + " " + 
				   std::to_string(radius*cos(i)*cos(alfa)) + "\n";
//			glVertex3f(radius*cos(i)*sin(alfa), radius*sin(i), radius*cos(i)*cos(alfa));
			
			ret += std::to_string(radius*cos(delta_a)*sin(alfa)) + " " + 
				   std::to_string(radius*sin(delta_a)) + " " + 
				   std::to_string(radius*cos(delta_a)*cos(alfa)) + "\n";
//			glVertex3f(radius*cos(delta_b)*sin(alfa), radius*sin(delta_b), radius*cos(delta_b)*cos(alfa));
			
			ret += std::to_string(radius*cos(delta_a)*sin(delta_a)) + " " + 
				   std::to_string(radius*sin(delta_b)) + " " + 
				   std::to_string(radius*cos(delta_b)*cos(delta_a)) + "\n";
//			glVertex3f(radius*cos(delta_b)*sin(delta_a), radius*sin(delta_b), radius*cos(delta_b)*cos(delta_a));
			
			ret += std::to_string(radius*cos(i)*sin(alfa)) + " " + 
				   std::to_string(radius*sin(i)) + " " + 
				   std::to_string(radius*cos(i)*cos(alfa)) + "\n";
//			glVertex3f(radius*cos(i)*sin(alfa), radius*sin(i), radius*cos(i)*cos(alfa));
			
			ret += std::to_string(radius*cos(delta_b)*sin(delta_a)) + " " + 
				   std::to_string(radius*sin(delta_b)) + " " + 
				   std::to_string(radius*cos(delta_b)*cos(delta_a)) + "\n";
//			glVertex3f(radius*cos(delta_b)*sin(delta_a), radius*sin(delta_b), radius*cos(delta_b)*cos(delta_a));
			
			ret += std::to_string(radius*cos(i)*sin(delta_a)) + " " + 
				   std::to_string(radius*sin(i)) + " " + 
				   std::to_string(radius*cos(i)*cos(delta_a)) + "\n";
//			glVertex3f(radius*cos(i)*sin(delta_a), radius*sin(i), radius*cos(i)*cos(delta_a));
	
			count+=6;
		}

	}

	return std::to_string(count) + "\n" + ret; 
}
/*
 *	TODO:	1: Make a class for vertices so that its easier to 
 *			   manipulate them (maybe not that necessary in 
 *			   this phase but almost definitely will be)
 *			
 *			2: Calculate the vertices (should we write them as
 *			   as we calculate them or do we accumalate them 
 *			   in some structure!?)
 *			
 *			3: About the Plane, the assignment doesnt talk about
 *			   dimensions.
 *			
 *			4: Every graphical primitive, where are the position?
 *			   how are we going to calculate the vertices if we
 *			   dont know the "origin" we're considering to draw
 *			   the damn things 
 */

/**
 *	./generator plane <lsize> plane.3d								| 4
 *	./generator box <x_dim> <y_dim> <z_dim> <n_divs>? box.3d		| 6/7
 *	./generator sphere <radius> <slices> <stacks> sphere.3d			| 6
 *	./generator cone <b_radius> <height> <slices> <stacks> cone.3d	| 7
 *
 *	this follows the following format:
 *		1* - name of the graphical primitive (only letters)
 *		2* - a series of numbers (digits and maybe '.')
 *		3* - finally the file name (digits, letters, '.', ... )
 *
 *	if we make it this way, we can verify the input using regex
 */

int main(int argc, char* argv[])
{
	if(argc < 2)
	{
		cout << "Wrong number of arguments!" << endl;
		return -1;
	}

	if(!validate_input(argc, argv))
	{
		cout << "Invalid arguments!" << endl;
		return -1;
	}

	if(!strcmp(argv[1], "plane"))
	{
		float info[3] = {-1.0f, -1.0f, -1.0f};
		for(int i = 2; i < argc-1; i++)
		{
			if(!strcmp("-d", argv[i]))
				info[2] = atof(argv[++i]);
			else if(info[0] < 0)
				info[0] = atof(argv[i]);
			else
				info[1] = atof(argv[i]);
		}
		generate_plane(info[0], info[1], info[2], argv[argc-1]);
	}

	if(!strcmp(argv[1], "sphere") && argc == 6)
	{
		ofstream file;
		file.open(argv[5]);
		file << drawSphere(atof(argv[2]), atof(argv[3]), atof(argv[4]));
		file.close();
	}

	return 0;
}

bool validate_input(int argc, char* argv[])
{
	/*
	 *	nums:	0.1 | .12 | 12 | 123.1
	 *
	 *	Should make this so that we always write to .3d files?
	 */
	regex nums("[0-9]*\\.?[0-9]+");
	regex flag("-[a-zA-Z]");
	regex file("[-_0-9a-zA-Z]+[.][0-9a-zA-Z]+");
	regex graph_prim("[a-zA-Z]+");
	bool numbers = true;
	bool f;

	bool gp = regex_match(argv[1], graph_prim);
	if(!gp) return false;

	for(int i = 2; i < argc-1 && numbers; i++ )
		if(regex_match(argv[i++], flag))
			numbers = regex_match(argv[i], nums);
	
	if(!numbers)
		return false;
	
	f = regex_match(argv[argc-1], file);

	return f;
}

void generate_sphere(char* argv[])
{
	/*	There are 4 elements in the argv argument
	 *
	 *	argv[0] = radius	(Valid floats | we already checked its validity)
	 *	argv[1] = slices	( same goes to this )
	 *	argv[2] = stacks	( and this )
	 *	argv[3] = file to write to	( and this )
	 *
	 *	So we can just convert the values we need to convert.
	 */

	 float radius = atof(argv[0]);
	 float slices = atof(argv[1]);
	 float stacks = atof(argv[2]);

	/*
	 *	Now for the real work:
	 *
	 *	I think we have to sort of "cross" the slices with the stacks
	 *	and "iterate" over the vertices it forms. I think looking at
	 *	the glutWireSphere can really help with visualizing this
	 *
	 */


}

void generate_plane(float x_len, float z_len, float n_divs ,char* file)
{
	cout << x_len << endl;
	cout << z_len << endl;
	cout << n_divs << endl;
}
