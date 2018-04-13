#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main(int argc, char *argv[])
{
    if(argc != 3) {
        printf("Usage: ./compare [tessellations path] [shortest path]\n");
        exit(EXIT_FAILURE);
    }
    FILE *tess = fopen(argv[1], "r");
    FILE *path = fopen(argv[2], "r");
    if(tess == NULL || path == NULL) {
        printf("File(s) not found. Exiting Program. Good Day.\n");
        exit(EXIT_FAILURE);
    }
    char buf[1024];
    int i = 0;
    int j = 0;
    int p1 = -1;
    int p2 = -1;
    int p3 = -1;
    int p4 = -1;
    int found = 0;

    int path_size = 0;
    while(fgets(buf, 1024, path)) {
        path_size++;
    }
    int tess_size = 0;
    while(fgets(buf, 1024, tess)) {
        tess_size++;
    }

    printf("\n... comparing %s and %s\n", argv[1], argv[2]);

    fclose(path);
    path = fopen(argv[2], "r");
    for(i = 0; i < path_size; i++) {
        fscanf(path, "%d %d", &p1, &p2);
        fclose(tess);
        tess = fopen(argv[1], "r");
        found = 0;
        for(j = 0; j < tess_size; j++) {
            fscanf(tess, "%d %d", &p3, &p4);
            if(((p1 == p3) && (p2 == p4)) || ((p2 == p3) && (p1 == p4))) {
                found = 1;
                break;
            }
        }
        if(!found) {
            printf("... PATH NOT FOUND\n\n");
            break;
        }
    }
    if(found) {
        printf("... PATH FOUND\n\n");
    }
    fclose(path);
    fclose(tess);
    return 0;
}
