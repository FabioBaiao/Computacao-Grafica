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
		file << drawSphere(...);
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
