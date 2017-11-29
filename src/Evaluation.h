/*
 * Evaluation.h
 *
 *  Created on: May 9, 2016
 *      Author: llx
 */

#ifndef EVALUATION_H_
#define EVALUATION_H_
#include <vector>
#include <ext/hash_map>
#include <sstream>
#include <cstdio>       // popen, printf, snprintf
#include <sys/wait.h>   // WIFEXITED() WEXITSTATUS()
#include <errno.h>      // extern int errno;
#include <string>
#include <cstring>
#include <cstdlib>
#include <cassert>
#define MAX_SIZE 1024
#define COMMAND_MAX_SIZE 1024

class Cgenome
{
public:
	std::fstream::pos_type index; //start from 0
	int length;
	Cgenome() :
			length(0), index(0)
	{
	}
	;
};

class Clongread
{
public:
	std::fstream::pos_type index; //start from 0
	int length;
	Clongread() :
			length(0), index(0)
	{
	}
	;
};

namespace __gnu_cxx
{
struct str_hash
{
	size_t operator()(const std::string& str) const;
};
struct str_equal
{
	bool operator()(const std::string& s1, const std::string& s2) const;
};
}

enum errtype
{
	substitution, insertion, deletion,
};

struct errnode
{
	unsigned long long location;
	errtype err;
};

struct mapline
{
	std::string longreadname;
	std::string genomename;
	unsigned long long startoffset;
	unsigned long long endoffset; //not included
	std::vector<errnode> errnodes;
	mapline():longreadname("*"){};
};

struct mapgroup
{
	std::string filename;
	std::vector<mapline> maplines;
};

bool myexec(const char *cmd, std::stringstream &resvec) {
    resvec.clear();
    FILE *pp = popen(cmd, "r");
    if (!pp) {
        return false;
    }
    char tmp[1024];
    while (fgets(tmp, sizeof(tmp), pp) != NULL) {
        resvec << tmp;
    }
    std::fstream logfile;
    logfile.open("log.txt", std::ios::app);
    logfile<< resvec.rdbuf();
    pclose(pp);
    return true;
}

bool myexec2(const char *command, std::stringstream &resvec,std::string& errmsg) {
	 	assert(command);

	    char buffer[MAX_SIZE] = {'\0'};
	    std::string final_msg;

	    // the exit status of the command.
	    int rc = 0;

	    // I/O redirection.
	    char cmd[COMMAND_MAX_SIZE] = {'\0'};
	    std::stringstream ss;
	    std::string tid;
	    ss << omp_get_thread_num();
	    ss >> tid;
	    std::string errfilename = "errs"+ tid +".log" ;
	    snprintf(cmd, sizeof(cmd), "%s 2>>%s", command,errfilename.c_str());

	    FILE *fp = popen(cmd, "r");
	    if(NULL == fp)
	    {   // if fork(2) or pipe(2) calls fail, or cannot callocate memory.
	        // it does not set errno if memory allocation fails.
	        // if the underlying fork(2) or pipe(2) fails, errno is set
	        // appropriately.
	        // if the type arguments is invalid, and this condition is detected,
	        // errno is set to EINVAL.
	        snprintf(buffer, sizeof(buffer),
	                "popen failed. %s, with errno %d.\n", strerror(errno), errno);
	        final_msg = buffer;
	        errmsg = final_msg.c_str();
	        return false;
	    }

	    char result[MAX_SIZE] = {'\0'};
	    std::string child_result;
	    while(fgets(result, sizeof(result), fp) != NULL)
	    {
	    	resvec << result;
	    }

	    // waits for the associated process to terminate and returns
	    // the exit status of the command as returned by wait4(2).
	    rc = pclose(fp);
	    if(-1 == rc)
	    {   // return -1 if wait4(2) returns an error, or some other error is detected.
	        // if pclose cannot obtain the child status, errno is set to ECHILD.

	        if(ECHILD==errno) {
	            final_msg += "pclose cannot obtain the child status.\n";
	        }
	        else {
	            snprintf(buffer, sizeof(buffer),
	                    "Close file failed. %s, with errno %d.\n", strerror(errno), errno);
	            final_msg += buffer;
	        }

	        errmsg = final_msg.c_str();
	        return false;
	    }

	    // child process exit status.
	    int status_child = WEXITSTATUS(rc);

	    // the success message is here.
	    final_msg += child_result;
	    snprintf(buffer, sizeof(buffer),
	            "[%s]: command exit status [%d] and child process exit status [%d].\r\n",
	            command, rc, status_child);
	    final_msg += buffer;
	    errmsg = final_msg.c_str();

	    if(status_child==0) {
	        // child process exits SUCCESS.
	        return true;
	    }
	    else {
	        // child process exits FAILED.
	        return false;
	    }
}

class Cfilebuffer
{
private:
	char* buffer;
	std::string myfilename;
	long mysize;
	long expectsize;
	long startoffset;
	long hited;
	long nothited;
public:
	Cfilebuffer():hited(0),nothited(0){};
	Cfilebuffer(std::string filename,long size):hited(0),nothited(0){
		buffer = new char[size+1];
		buffer[size] = '\0';
		myfilename = filename;
		expectsize = size;
		startoffset = 0;
		std::ifstream file(filename.c_str());
		file.read(buffer,expectsize);
		mysize = file.gcount();
		buffer[mysize]='\0';
		file.close();
	}
	void refreshbuffer(std::string filename,long size){
			delete[] buffer;
			buffer = new char[size+1];
			buffer[size] = '\0';
			myfilename = filename;
			expectsize = size;
			startoffset = 0;
			std::ifstream file(filename.c_str());
			file.read(buffer,expectsize);
			mysize = file.gcount();
			buffer[mysize]='\0';
			file.close();
	}
	std::string Getstring(char *argv, std::fstream::pos_type begin, std::fstream::pos_type end);
	~Cfilebuffer()
	{
		delete[] buffer;
	}
	float hitraio()
	{
		return (float)hited/(hited+nothited);
	}
};

#endif /* EVALUATION_H_ */
