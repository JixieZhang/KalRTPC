/*Class to read and store the configuration*/

#include "ConfigManager.hh"


ConfigManager* ConfigManager::fConfigManager=0;
ConfigManager* ConfigManager::GetConfigManager()
{ 
	if(!fConfigManager) 
	{
		std::cout<<"No instance of ConfigManager found! "
			<<"Use ConfigManager::ConfigManager(const char* filename) to build one!"<<std::endl;
	}
	return fConfigManager; 
}
//..............................................................................
ConfigManager::ConfigManager(const char* filename)
{
#ifdef USAGEMANAGER_DEBUG
	if(Global_Debug_Level < (int)USAGEMANAGER_DEBUG) 
	{
		Global_Debug_Level = (int)USAGEMANAGER_DEBUG;
		std::cout<<"set Global_Debug_Level to "<<Global_Debug_Level<<std::endl;
	}
#endif

	if(fConfigManager)
	{ 
		cout<<"Error! fConfigManager constructed twice.\n"
			<<"Use ConfigManager::GetConfigManager() to get the static pointer \n"; 
	}
	fConfigManager=this;

    bConfigLoaded=false;
    bConfigLoaded=this->ReadFile(filename);
}

//..............................................................................
ConfigManager::~ConfigManager()
{    
	MapParam_s.clear();

	cout << "deleting fConfigManager ......" << endl;
	fConfigManager = 0;
}


//////////////////////////////////////////////////////////////////////
//remove space or tab \t at both end of the string
void  ConfigManager::Trim(std::string &astr)
{
	if(astr.length()<1) return;
	size_t i=astr.length()-1;
	while(astr[i]==' ' || astr[i]=='\t' || astr[i]=='\r' )
	{
		astr.erase(astr.length()-1,1);
		i=astr.length()-1;
	}
	i=0;
	while(astr[i]==' ' || astr[i]=='\t' || astr[i]=='\r' )
	{
		astr.erase(0,1);
	}
}

void ConfigManager::Trim(char *astr)
{
	if(strlen(astr)<1) return;

	size_t i=strlen(astr)-1;
	bool search=true;
	do{
		if(astr[i]==' ' || astr[i]=='\t' || astr[i]=='\r' ) 
		{
			astr[i]='\0';
			i--;
		}
		else search=false;
	}while(search);

	i=0;
	search=true;
	do{
		if(astr[i]==' ' || astr[i]=='\t' || astr[i]=='\r' ) 
		{
			i++;
		}
		else search=false;
	}while(search);

	if(i>0)
	{
		std::string cppstr(astr+i);
		strcpy(astr,cppstr.c_str());
	}
}

void ConfigManager::PrintParamMap()
{
	if(MapParam_s.size() <= 0) return;
	cout<<"*************************MAP**VALUE**START***************************"<<endl;
	cout<<"PrintParamMap(): The content of the parameters in the map are:\n";
	map<string,string>::iterator it_s;
	for ( it_s=MapParam_s.begin() ; it_s != MapParam_s.end(); it_s++ )
	{
		cout<<"MapParam_s["<<setw(12)<<(*it_s).first<<"] = \""<<(*it_s).second<<"\""<<endl;
	}
	cout<<"*************************MAP**VALUE****END***************************"<<endl;
}

//c++ version of function to read a parameter file
//Input: filename c string
//Output: return false if file can not open succesfully, otherwise return true 
//        all parameters written in the given file will be stored and get be
//        achieved by calling bool GetParameter<typename>(char *name, T &value);
//        Please NOTE that the typename could be size_t,int,float,double,bool,char*,string only
//feature: 
//1)The parameter files must in a format as "parameter=value", allow space ' ' beside '='. 
//2)The value is ended by comma ',', semi-comma ';' or space ' ', care should be taken for string parameters
//2)Use '#' to comment out lines or part of lines. Everything after # will be ignored.
//3)A line without '=' will not be understood. All non-understood lines will be printed out and ignored.
//4)The maximum length of each line is 1024.

bool ConfigManager::ReadFile(const char *filename)
{
	ifstream ini(filename,ios_base::in);
	if(!ini.good())
	{
		cout<<"***Error! Can not open parameter file \""<<filename<<"\"!"<<endl;
		return false;
	}
	else  cout<<"Reading parameter file \""<<filename<<"\""<<endl;

	int nline=0;
	size_t found=0;
	static const int MaxLineSize=1024;
	char line[MaxLineSize];
	std::string strLine;
	std::string strName;
	std::string strValue;

	while(!ini.eof())
	{
		ini.getline(line,MaxLineSize);
		nline++;
		//Search forward. for the '#' to skip the afterwards
		strLine=line;
		found=strLine.find('#');
		if(found!=string::npos)	strLine=strLine.substr(0,found);
	
		//trim the string
		Trim(strLine);
		if(strLine.length()==0) continue; //this is a comment line whose first char is '#'

		//try to find '=', if not find, ignore this line		
		found=strLine.find('=');
		if (found==string::npos || strLine.length()<3) 
		{
			cout<<"Warning! line "<<nline<<" is not recognized: \""<<line<<"\" \n";
			continue;
		}
		strName=strLine.substr(0,found);
		Trim(strName);	

		strValue=strLine.substr(found+1);
		found=strValue.find(';');
		if (found!=string::npos) strValue=strValue.substr(0,found);
		Trim(strValue);	
		cout<<setw(15)<<strName<<" = "<<strValue<<endl;
		//store the value into the map
		MapParam_s[strName]=strValue;

	}
	ini.close();

	return true;

}

bool ConfigManager::ReadFile_c(const char *filename)
{
	FILE *ini;
	if((ini=fopen(filename,"r"))==NULL)
	{
		printf("***Error! Can not open parameter file \"%s\"!\n",filename);
		return false;
	}
	else  printf("Reading parameter file \"%s\"\n",filename);

	char ch[]="=;\n";
	int nline=0;
	static const int MaxLineSize=1024;
	char *name,*value,line[MaxLineSize];
	char *pDest;
	std::string strName,strValue;

	while(!feof(ini))
	{
		//empty the buffer
		sprintf(line,"%s","");
		fgets(line,MaxLineSize,ini);
		nline++;

		//Search forward. for the '#' to skip the fterwards
		pDest = strchr( line, '#' );
		if(pDest!=NULL) *pDest='\0'; 
		//trim the string
		Trim(line);
		if(strlen(line)==0) continue; //this is a comment line

		//try to find '=', if not find, ignore this line		
		pDest = strchr( line, '=' );
		if(pDest==NULL || strlen(line)<3) 
		{
			printf("Warning! line %d is not recognized: \"%s\" \n",nline,line);
			continue;
		}

		name=strtok(line,ch);
		value=strtok(0,ch);

		//store the value into the map
		Trim(name);
		Trim(value);
		strName=name;
		strValue=value;
		printf("%15s = %s\n",strName.c_str(),strValue.c_str());
		MapParam_s[strName]=strValue;
		
	}
	fclose(ini);

	return true;
}


//using the branch vector to find the variables
std::string ConfigManager::GetParameter(std::string name)
{	
	bool found=false;
	std::string parameter("");
	map<string,string>::iterator it;

	/*
	//brute force method to find the key
	for ( it=MapParam_s.begin() ; it != MapParam_s.end(); it++ )
	{
	if((*it).first.compare(name)==0) 
	{
	parameter=(*it).second; 
	found=true;
	break;
	}
	}
	*/
	it=MapParam_s.find(name);
	if(it != MapParam_s.end()) 
	{
		parameter=(*it).second; 
		found=true;
	}

	if(!found) std::cout<<"GetParameter(Name="<<name<<"): wrong parameter name!!! "<<std::endl;
	return parameter;
}

std::string ConfigManager::GetParameter(const char *name)
{ 	
	std::string para(name);
	return GetParameter(para);
}

/*
bool ConfigManager::GetParameter(const std::string name,size_t &parameter)
{
std::string para=GetParameter(name);
bool found=(para.length()==0)?false:true;
parameter=(size_t) atoi(para.c_str());
return found;
}
bool ConfigManager::GetParameter(const std::string name,int    &parameter)
{
std::string para=GetParameter(name);
bool found=(para.length()==0)?false:true;
parameter=(int) atoi(para.c_str());
return found;
}
bool ConfigManager::GetParameter(const std::string name,float  &parameter)
{
std::string para=GetParameter(name);
bool found=(para.length()==0)?false:true;
parameter=(float) atof(para.c_str());
return found;
}
bool ConfigManager::GetParameter(const std::string name,double &parameter)
{
std::string para=GetParameter(name);
bool found=(para.length()==0)?false:true;
parameter=(double) atof(para.c_str());
return found;
}

*/
bool ConfigManager::GetParameter(const std::string name,std::string &parameter)
{
	parameter=GetParameter(name);
	bool found=(parameter.length()==0)?false:true;
	return found;
}
bool ConfigManager::GetParameter(const std::string name,char *parameter)
{
	std::string para=GetParameter(name);
	bool found=(para.length()==0)?false:true;
	sprintf(parameter,"%s",para.c_str());
	return found;
}

//////////////////////////////////////////////////////////////////

bool ConfigManager::SetParameter(std::string name,std::string value)
{
	bool found=false;
	std::string parameter("");
	map<string,string>::iterator it;

	it=MapParam_s.find(name);
	if(it != MapParam_s.end()) 
	{
		parameter=(*it).second; 
		found=true;
		std::cout<<"###SetParameter(name="<<name<<", valee="<<value
			<<")  called!###"<<std::endl;
		it->second=value; 
	}

	if(!found) 
	{
		std::cout<<"SetParameter(name="<<name<<", valee="<<value
			<<"): wrong parameter name!!! exit ... "<<std::endl;
		exit(-999);
	}
	return found;
}