#define TIXML_USE_STL 
#define BUFF_SIZE 128
#include "tinyxml/tinyxml.h"
#include "point.h"
#include <vector>
#include <fstream>
#include <string>

using namespace std;

//We assume that the .xml and .3d files passed are correct.
string parseFileName(string filenameFromXml){
	// for some reason, the '\"' symbol uses 3bytes
	return filenameFromXml.substr(3, filenameFromXml.length()-6); 
}

int main(int argc, char** argv){
	if(argc != 2){
		cout << "The command invocation must be in the form 'engine filename'\n";
		return -1;
	}

	TiXmlDocument doc(argv[1]);
	bool loadOkay = doc.LoadFile();
	if(!loadOkay){
		cout << "Error to load file " << argv[1] << ".\n";
		return -1;
	}

	vector<point> points{};
	TiXmlHandle docHandle(&doc);
	TiXmlElement* model = doc.FirstChild("scene")->FirstChild("model")->ToElement();
	for(; model; model=model->NextSiblingElement()){
		string filename; 
		int r = model->QueryStringAttribute("file", &filename);
		if(r == TIXML_SUCCESS){
			parseFileName(filename);
			ifstream file;
			file.open(filename);
			int n_vertex; 
			file >> n_vertex;
			cout << n_vertex;
		}
	}
}
