/*Class to read and store the configuration*/
#ifndef _CONFIGMANAGER_H_
#define _CONFIGMANAGER_H_

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <string>
#include <vector>
#include <map>
#include <typeinfo>

#include <ctype.h>
using namespace std;

class ConfigManager
{
public: 
	//Static method which returns the singleton pointer of ConfigManager or its derived class.
	static ConfigManager* GetConfigManager();

private:
	static ConfigManager* fConfigManager;

public:
    ConfigManager(const char* filename="Detector.ini");
    ~ConfigManager();

	//print the parameter map, which is the name and values written in the parameter files
	void PrintParamMap();

	bool SetParameter(std::string name,std::string value);

	//The majoy GET function for the parameters which were written in the parameter file
	//Get the parameter by name, which can be found in the parameter file, return an empty string if not found
	std::string GetParameter(std::string name);
	std::string GetParameter(const char *name);
	
	//Get argument which has the same type as the given one, return true if found otherwise return false 
	bool GetParameter(const std::string name,std::string &argument);
	bool GetParameter(const std::string name,char *argument);
	//the following five have been included in the template
	//bool GetParameter(const std::string name,size_t &argument);
	//bool GetParameter(const std::string name,int    &argument);
	//bool GetParameter(const std::string name,float  &argument);
	//bool GetParameter(const std::string name,double &argument);
	//bool GetParameter(const std::string name,bool  &argument);
	
	//Get argument which has the same type as the given one, return true if found otherwise return false 
	//Basiclly all regular 'value' type, i.e signed and unsigned short, int, long, float and double
	//will be included in this template
	template <typename T>
	bool GetParameter(const std::string name,T &argument)
	{
		bool found=true;
		std::string arg=GetParameter(name);
		if(arg.length()<1) return false;

		if(typeid(argument)==typeid(float)    || typeid(argument)==typeid(double)) 
			argument=(T)atof(arg.c_str());
		else if(typeid(argument)==typeid(int) || typeid(argument)==typeid(unsigned int) ||
			typeid(argument)==typeid(long)    || typeid(argument)==typeid(unsigned long) ||
			typeid(argument)==typeid(short)   || typeid(argument)==typeid(unsigned short) ) 
			argument=(T)atoi(arg.c_str());
		else if(typeid(argument)==typeid(bool)) 
			argument=(arg.compare("true")==0  || arg.compare("TRUE")==0) ? true : false;
		else 
		{
			std::cout<<"unknown input type \""<<typeid(argument).name()<<std::endl;
			found=false;
		}
		return found;
	}


	//c++ version of function to read a parameter file
	//Input: filename c string
	//Output: return false if file can not open succesfully, otherwise return true 
	//        all parameters written in the given file will be stored and get be
	//        achieved by calling bool GetParameter<typename>(char *name, T &value);
	//        Please NOTE that the typename could be size_t,int,float,double,bool,char*,string only
	//feature: 
	//1)The parameter files must in a format as "parameter=value", allow space ' ' beside '='. 
	//2)The value is ended by semi-comma ';'. Care should be taken for string parameters
	//2)Use '#' to comment out lines or part of lines. Everything after # will be ignored.
	//3)A line without '=' will not be understood. All non-understood lines will be printed out and ignored.
	//4)The maximum length of each line is 1024.
	bool ReadFile_c(const char *filename);

	//c version of function to read a parameter file
	//Input: filename c string
	//Output: return false if file can not open succesfully, otherwise return true 
	//        all parameters written in the given file will be stored and get be
	//        achieved by calling bool GetParameter<typename>(char *name, T &value);
	//        Please NOTE that the typename could be size_t,int,float,double,bool,char*,string only
	//feature: 
	//1)The parameter files must in a format as "parameter=value", allow space ' ' beside '='. 
	//2)The value is ended by semi-comma ';'. Care should be taken for string parameters
	//2)Use '#' to comment out lines or part of lines. Everything after # will be ignored.
	//3)A line without '=' will not be understood. All non-understood lines will be printed out and ignored.
	//4)The maximum length of each line is 1024.
	bool ReadFile(const char *filename);

private:
	//remove space or tab \t at both ends of the string
	static void Trim(std::string &astr);
	static void Trim(char *astr);

	bool bConfigLoaded;
	map<string,string> MapParam_s;

};

typedef ConfigManager ConfigT;
typedef ConfigManager ReadConfig;

#endif //_CONFIGMANAGER_H_
