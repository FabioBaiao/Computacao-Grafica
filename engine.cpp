#define TIXML_USE_STL 
#include "tinyxml/tinystr.h"
#include "tinyxml/tinyxml.h"

int main(int argc, char** argv){
	if(argc != 2){
		std::cout << "The command invocation must be in the form 'engine filename'\n";
		return -1;
	}

	TiXmlDocument doc(argv[1]);
	bool loadOkay = doc.LoadFile();
	if(!loadOkay){
		std::cout << "Error to load file " << argv[1] << ".\n";
		return -1;
	}
}
