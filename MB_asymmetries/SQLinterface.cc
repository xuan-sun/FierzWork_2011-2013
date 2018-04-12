#include "SQLinterface.hh"
#include <stdlib.h>
#include <stdio.h>

SQLdatabase::SQLdatabase(const std::string& name,
			 const std::string& dbAddress,
			 const std::string& dbUser,
			 const std::string& dbPass,
			 unsigned int port,
			 unsigned int ntries): numFields(0), db(NULL),res(NULL),dbname(name) {

  char temp[10];
  sprintf(temp,"%i",port);
  std::string dbAddressFull = "mysql://"+dbAddress+":"+std::string(temp)+"/"+dbname;
  
  while (!db)
    {
      ntries--;
      db = TSQLServer::Connect(dbAddressFull.c_str(),dbUser.c_str(),dbPass.c_str());
      if(!db) {
	if(!ntries)
	  break;
	printf("** DB Connection %s@%s failed... retrying...\n",dbUser.c_str(),dbAddressFull.c_str());
	//sleep(2);
      } else {
	printf("Connected to DB server: %s\n", db->ServerInfo());
	return;
      }
    }

  if (!db) throw "Couldn't Connect to Database";
}
  
void SQLdatabase::fetchQuery(const char* q) {
  if (!db) {res=NULL;}
  else {
    if (!q) q=query.c_str();
    if (res) delete(res);
    res = db->Query(q);
    
    if(db->GetErrorCode()) {
      throw "DB query Error. Check Query Syntax";
    }
    else numFields = res->GetFieldCount();
  }
}

bool SQLdatabase::queryReturnsTrue() {
  if (!res) return false;
  else if (res->GetFieldCount()>0) return true;
  else return false;
};
    

std::string SQLdatabase::returnQueryEntry(int* field, TSQLResult* r) {
  std::string result;
  if (!r) r=res;
  else numFields = r->GetFieldCount();
  TSQLRow *row = r->Next();
  if (!field) {
    if (!row) return "";// "Row doesn't exist in Query. Query may be empty";
    for (int f = 0; f<numFields; f++)
      {
	result += std::string(row->GetField(f));
	if (f!=(numFields-1)) result += " ";
      }
  }
  else {
    if (*field>=numFields) throw "Chose field outside number of fields in query";
    result += std::string(row->GetField(*field));
  }
  delete row;
  return result;
}
  

