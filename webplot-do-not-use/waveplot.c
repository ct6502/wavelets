#include <stdio.h>
#ifndef NO_STDLIB_H
#include <stdlib.h>
#else
char *getenv();
#endif
#include <string.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <time.h>

#define MAX_ENTRIES 2500

typedef struct {
	char *name;
	char *val;
} entry;

char *makeword(char *line, char stop);
char *fmakeword(FILE *f, char stop, int *len);
char x2c(char *what);
void unescape_url(char *url);
void plustospace(char *str);

FILE  *infp,*outfp;

main(int argc, char *argv[]) {
	entry entries[MAX_ENTRIES];
	register int x,m=0;
	int cl;
	char data[255],filename[255],file2[255],dir[255],command[255],time_str[20];
	clock_t clocktime;
	time_t curtime;
	struct tm *loctime;
	
	printf("Content-type: text/html%c%c",10,10);

	curtime = time (NULL); /* Get the current time. */
	loctime = localtime (&curtime); /* Convert to localtime representation.*/
	strftime (time_str,20,"%H%M%S",loctime);  /* Convert to string */

/* local directories */
	sprintf(dir, "/export/home/httpd/htdocs/research/wavelets/plot/");
	sprintf(filename, "out/wave_%s", getenv("REMOTE_HOST"));

	if(strcmp(getenv("REQUEST_METHOD"),"POST")) {
		goto skip;
/*		printf("This script should be referenced using ");
		printf("<A HREF=\"waveplot.html\">Wavelet Plot</A>.%c",10);
		exit(1);
*/
	}

/* not sure what this does... */
	if(strcmp(getenv("CONTENT_TYPE"),"application/x-www-form-urlencoded")) {
		printf("<BODY>This script can only be used to decode form results. \n");
		exit(1);
	}
	cl = atoi(getenv("CONTENT_LENGTH"));

/* decode the form data */
	for(x=0;cl && (!feof(stdin));x++) {
		m=x;
		entries[x].val = fmakeword(stdin,'&',&cl);
		plustospace(entries[x].val);
		unescape_url(entries[x].val);
		entries[x].name = makeword(entries[x].val,'=');
	}

/* read in options */
	if (strcmp(entries[1].name,"$dataset")==0) {
		sprintf(command,"rm -f %s%s.opt",dir,filename);
		system(command);
		sprintf(file2,"%s%s.dat",dir,filename);}
	else
		sprintf(file2,"%s%s.opt",dir,filename);

/* remove previous gif & postscript files */
	sprintf(command,"rm -f %s%s*.gif %s%s*.ps",dir,filename,dir,filename);
	system(command);

/* write out options */
	outfp = fopen(file2,"w");
	for(x=0; x <= m; x++)
		fprintf(outfp,"%s\n%s\n",entries[x].name,entries[x].val);
	fclose(outfp);
	
skip:
/* construct IDL batch file */
	sprintf(file2,"%s%s.idl",dir,filename);
	outfp = fopen(file2,"w");
	fprintf(outfp,"CATCH,error_status\n");
	fprintf(outfp,"IF (error_status NE 0) THEN EXIT\n");
	fprintf(outfp,"CD,'%s'\n",dir);
	fprintf(outfp,"RESTORE,'waveplot.sav'\n");
	fprintf(outfp,".run waveplot\n");
	fprintf(outfp,"waveplot,'%s','%s'\n",filename,time_str);
	fprintf(outfp,"EXIT\n");
	fclose(outfp);

/* run IDL batch file */
	sprintf(command,"/usr/local/rsi/idl/bin/idl %s",file2);
	system(command);

/* read in HTML & dump to browser */
	sprintf(file2,"%s%s.html",dir,filename);
	infp = fopen(file2,"r");
	(void) fgets(data,sizeof(data),infp);
	while (!feof(infp)) {
		printf("%s\n",data);
		(void) fgets(data,sizeof(data),infp);
	}
	fclose(infp);

/* automatically remove data & gif files more than 4 days old... */
	sprintf(command,"find %sout -name 'wave_*' -mtime +4 -exec rm {} \\;",dir);
	system(command);

}
