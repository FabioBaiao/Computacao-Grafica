#define TIXML_USE_STL 
#include "tinyxml/tinyxml.h"
#include "point.h"
#include <vector>
#include <string>

//We assume that the .xml and .3d files passed are correct.
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

	std::vector<point> points{};
	TiXmlHandle docHandle(&doc);
	TiXmlElement* model = doc.FirstChild("scene")->FirstChild("model")->ToElement();
	for(; model; model=model->NextSiblingElement()){
		std::string filenameFromXml, filename; 
		int r = model->QueryStringAttribute("file", &filenameFromXml);
		if(r == TIXML_SUCCESS){
			filename = filenameFromXml.substr(1, filenameFromXml.length()-2);
			//ifstream f = open(filename); 
			std::cout << filename;
		}
	}
}
