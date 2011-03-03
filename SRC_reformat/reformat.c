#include <stdio.h>

main(argc,argv)
  int argc;
  char *argv[];
{ double xxx=1;
  int i=0, zaehler,n;
  int x;
  double a,b;

 /* if (argc!=2) {
    puts("Ein Argument: Anzahl der Zahlen vor Newline");
    exit(-1); }*/

  scanf("%d%d%le%le",&x, &zaehler,&a,&b);
  
  /*zaehler=atoi(argv[1]);*/

  n=scanf("%le",&xxx);
  while(n!=-1){
    printf("%lg\n",xxx);
    i++;
    if (i >= zaehler) { 
      i=0;
      puts(""); }
    n=scanf("%lg",&xxx); }
}
