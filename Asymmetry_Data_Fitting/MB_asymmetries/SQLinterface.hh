#ifndef SQL_UTILS_HH
#define SQL_UTILS_HH


#include <TSQLServer.h>
#include <TSQLRow.h>
#include <TSQLResult.h>
#include <vector>
#include <string>

class SQLdatabase {
public:
  //Constructor
  SQLdatabase(const std::string& name,
	      const std::string& dbAddress,
	      const std::string& dbUser,
	      const std::string& dbPass,
	      unsigned int port = 3306,
	      unsigned int ntries = 3);

  //Destructor
  ~SQLdatabase() {if (res) delete (res); if (db) db->Close();}

  void fetchQuery(const char* q=NULL); //run an SQL query and assign result to res
  
  bool queryReturnsTrue();

  std::string returnQueryEntry(int* field = NULL, TSQLResult* r = NULL); //This will return a string with each field value separated by a space, unless a particular field number is declared for int field, in which case it returns that field value only. You will need to convert this string to whatever type suits your use

  int numFields;      // number of columns in a query
  std::string query;	//can be set to some default query if desired

protected:
  TSQLServer *db; //database being used
  TSQLResult *res; //Result of a query
  std::string dbname;


};

#endif
