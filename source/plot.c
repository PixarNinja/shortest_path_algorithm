#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[])
{
    if(argc != 2) {
        printf("Usage: ./plot [path to gnu_plot commands]\n");
        exit(EXIT_FAILURE);
    }
    char *plot = malloc(strlen("gnuplot -persistent ") + strlen(argv[1]) + 1);
    sprintf(plot, "gnuplot -persistent %s", argv[1]);
    system(plot);
    return 0;
}
