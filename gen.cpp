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

int main(int argc, char ** argv){
	if(argc == 1){
		print_help();
		return 0;
	} 

	std::string primitive(argv[1]); 
	std::vector<Point> points;
	char* filename;

	if(primitive == "plane" && argc == 5){
		//validate arguments
		ofstream file = open(argv[4]);
		file << drawPlane(...);
		close(file);
	}
	else if(primitive == "box" && argc == 7){
		//validate arguments
		ofstream file = open(argv[6]);
		file << drawBox(...);
		close(file);
	}
	else if(primitive == "sphere" && argc == 6){
		//validate arguments
		ofstream file = open(argv[5]);
		file << drawSphere(atof(argv[2]), atof(argv[3]), atof(argv[4]));
		close(file);
	}
	else if(primitive == "cone" && argc == 7){
		//validate arguments
		ofstream file = open(argv[6]);
		file << drawCone(...);
		close(file);
	}
	else if(primitive == "torus" && argc == 10000){
		//to-do if time available
	}
	else if(primitive == "cilinder" && argc == 1000){
		//to-do if time available
	}
	else{
		print_help();
		return -1;
	}
}
