#include <stdio.h>
#include <string.h>
#include <stdlib.h>
// #include <my_global.h>
#include </pub/include/mysql/mysql.h>
// #include <my_getopt.h>
#include <sys/syscall.h>
#include "variant.h"

#define MyBufLen 10240
char myBuffer[MyBufLen + 1]; 

MYSQL * mysqlConnection, mysqlMysql;
MYSQL_RES * mysqlResult, * mysqlResult0;
MYSQL_ROW mysqlRow, mysqlRow0;

int executeSQL(const char * sql)
{
   int state;
   /* create tables */
   /*state = mysql_query(mysqlConnection, sql);*/
   state = mysql_query(mysqlConnection,sql);
   if(state)
   {
      printf("MySQL error[%d]: %s\n", state, mysql_error(mysqlConnection));
      return 0;
   }
   return 1;
}

/*****************************************/
int connectDB(char * hostName, char * dbName, char * user, char * passwd)
{
   int state;
   /******************************************************
    *     * connect to MySQL
    *         *****************************************************/
   mysql_init(&mysqlMysql);
   mysqlConnection = mysql_real_connect(&mysqlMysql, hostName, user, passwd, dbName, 0, "/var/lib/mysql/mysql.sock", 0);
   if(mysqlConnection == NULL)
   {
      printf("DB %s connection error: %s\n>", dbName, mysql_error(&mysqlMysql));
      return 0;
   }
   else
      {
          printf(">>> Connected to the DB %s\n", dbName);
          /* use the DB */
          sprintf(myBuffer, "USE %s", dbName);
          state = mysql_query(mysqlConnection, myBuffer);
          if(state)
          {
             printf("    DB: %s cannot be used.\n", dbName);
             return 0;
          }
          else
          printf("   Changed to DB: %s.\n", dbName);
      }
    return 1;
}

void disconnectDB(const char * dbName)
{
   /* disconnected from MySQL */
   mysql_close(mysqlConnection);
   printf("<<< Disconnected from the DB %s\n", dbName);
}

void pingAndReconnectToDB()
{
   while(mysql_ping(&mysqlMysql))
   printf("   MySQL error: attempting reconnection...\n");
}

//////////////////////////////////////////////////////////////////////////////

