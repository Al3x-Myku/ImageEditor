#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <ctype.h>
#define true 1
#define false 0
#define ALLOCATE_MATRIX(matrix, rand, size) \
    (matrix) = (double **)malloc((rand) * (size))
#define ALLOCATE_ROW(matrix, coloane, size) \
    (matrix) = (double *)malloc((coloane) * (size))
#define FREE_2D_ARRAY(arr, rand) do { \
    for (int i = 0; i < (int)rand; i++) { \
        free(arr[i]); \
    } \
    free(arr); \
} while(0)


#define BASIC_SIZE 100
#define S1 sizeof(double)
#define S2 sizeof(double *)

typedef struct check {
	uint32_t e_loaded;
	uint32_t e_all;
	uint32_t e_rgb;
	uint32_t e_first;
} check;

typedef struct image {
	char *magic_nr;
	check state;
	uint32_t first_selection;
	uint32_t selection_counter;
	uint32_t *coord;
	uint32_t coloane;
	uint32_t rand;
	uint32_t intensity;
	double **img;
	double **mat_r;
	double **mat_g;
	double **mat_b;
} image;

typedef struct nucleu {
	int *edge;
	int *blur;
	int *gauss;
	int *sharpen;
} nucleu;

void initializare(image *ph)
{
    *ph = (image){0};
}

void initializare_filt(nucleu *filt)
{
 
    int *filters[4];
    for (int i = 0; i < 4; ++i) {
        filters[i] = (int *)malloc(9 * sizeof(int));
    }
    filt->edge = filters[0];
    filt->blur = filters[1];
    filt->gauss = filters[2];
    filt->sharpen = filters[3];

   
    int edge[9] = {
        -1, -1, -1,
        -1,  8, -1,
        -1, -1, -1
    };
    int sharpen[9] = {
         0, -1,  0,
        -1,  5, -1,
         0, -1,  0
    };
    int blur[9] = {
         1,  1,  1,
         1,  1,  1,
         1,  1,  1
    };
    int gauss[9] = {
         1,  2,  1,
         2,  4,  2,
         1,  2,  1
    };

   
    memcpy(filt->edge, edge, 9 * sizeof(int));
    memcpy(filt->sharpen, sharpen, 9 * sizeof(int));
    memcpy(filt->blur, blur, 9 * sizeof(int));
    memcpy(filt->gauss, gauss, 9 * sizeof(int));
}

void loading(image **ph, int type1, int type2, int pos, char *file)
{
    FILE *f = fopen(file, type1 ? "rt" : "rb");
    if (!f) return; 

    fseek(f, pos + 1, SEEK_SET);

    int rand = (*ph)->rand;
    int coloane = (*ph)->coloane;

    if (type2) {
        ALLOCATE_MATRIX((*ph)->img,rand, S2);
        for (int i = 0; i < rand; i++) {
            ALLOCATE_ROW((*ph)->img[i],coloane,S1);
            for (int j = 0; j < coloane; j++) {
                uint8_t tmp_ch;

                if (type1) {
                    fscanf(f, "%hhu", &tmp_ch);
                } else {
                    fread(&tmp_ch, sizeof(char), 1, f);
                }

                int tmp = (int)tmp_ch;
                tmp = (tmp > 255) ? tmp % 255 - 1 : (tmp < 0) ? 255 + tmp % 255 + 1 : tmp;

                (*ph)->img[i][j] = (double)tmp;
            }
        }
        (*ph)->state.e_rgb = false;
    } else {
        	ALLOCATE_MATRIX((*ph)->mat_r, rand, S2);
			ALLOCATE_MATRIX((*ph)->mat_g, rand, S2);
			ALLOCATE_MATRIX((*ph)->mat_b, rand, S2);


        for (int i = 0; i < rand; i++) {
    		ALLOCATE_ROW((*ph)->mat_r[i], coloane, S1);
  		 	 ALLOCATE_ROW((*ph)->mat_g[i], coloane, S1);
   			 ALLOCATE_ROW((*ph)->mat_b[i], coloane, S1);

            for (int j = 0; j < coloane; j++) {
                unsigned char tmp_ch[3];

                if (type1) {
                    fscanf(f, "%hhu %hhu %hhu", &tmp_ch[0], &tmp_ch[1], &tmp_ch[2]);
                } else {
                    fread(tmp_ch, sizeof(unsigned char), 3, f);
                }

                for (int k = 0; k < 3; k++) {
                    int tmp = (int)tmp_ch[k];
                    tmp = (tmp > 255) ? tmp % 255 - 1 : (tmp < 0) ? 255 + tmp % 255 + 1 : tmp;

                    if (k == 0) (*ph)->mat_r[i][j] = (double)tmp;
                    else if (k == 1) (*ph)->mat_g[i][j] = (double)tmp;
                    else (*ph)->mat_b[i][j] = (double)tmp;
                }
            }
        }
        (*ph)->state.e_rgb = true;
    }

    (*ph)->state.e_all = true;
    (*ph)->state.e_loaded = true;
    fclose(f);
}

void reset(image **ph) {
    image *img = *ph;

    if (img->state.e_rgb) {
        FREE_2D_ARRAY(img->mat_r, img->rand);
        FREE_2D_ARRAY(img->mat_g, img->rand);
        FREE_2D_ARRAY(img->mat_b, img->rand);
    } else { 
        FREE_2D_ARRAY(img->img, img->rand);
    }

    if (img->magic_nr) {
        free(img->magic_nr);
    }

    if (img->state.e_first) {
        free(img->coord);
    }

    *img = (image){
		.state.e_all = 0,
		.state.e_first = 0,
		.state.e_loaded = 0,
		.state.e_rgb = 0,
        .selection_counter = 0,
        .first_selection = 0
    };
}

void incaraca(image *ph, char *reading_file)
{   FILE *f = fopen(reading_file, "rt");
		if (ph->state.e_loaded || !f) {
			if (!f) {
				if (ph->state.e_loaded)
					reset(&*&ph);
				printf("Failed to load %s\n", reading_file);
				return;
			}
			if (strcmp(ph->magic_nr, "P2") == 0 || strcmp(ph->magic_nr, "P5") == 0) {
    FREE_2D_ARRAY(ph->img, ph->rand);
}

if (strcmp(ph->magic_nr, "P3") == 0 || strcmp(ph->magic_nr, "P6") == 0) {
        FREE_2D_ARRAY(ph->mat_r, ph->rand);
        FREE_2D_ARRAY(ph->mat_g, ph->rand);
        FREE_2D_ARRAY(ph->mat_b, ph->rand);
}

free(ph->magic_nr);

		}

	ph->magic_nr =malloc(3);
	char *m_ptr = ph->magic_nr;
fgets(m_ptr, 3, f);

int *cols = (int*)&(ph->coloane), *rows = (int*)&(ph->rand), *inten = (int*)&(ph->intensity);
fscanf(f, "%d%d", cols, rows);
fscanf(f, "\n%d", inten);

long current_pos = ftell(f);
int pos = (int)current_pos;

	fclose(f);
 if (strcmp(ph->magic_nr, "P2") == 0) {
        loading(&ph, 1, 1, pos, reading_file);
    } else if (strcmp(ph->magic_nr, "P3") == 0) {
        loading(&ph, 1, 0, pos, reading_file);
    } else if (strcmp(ph->magic_nr, "P5") == 0) {
        loading(&ph, 0, 1, pos, reading_file);
    } else if (strcmp(ph->magic_nr, "P6") == 0) {
        loading(&ph, 0, 0, pos, reading_file);
    } else {
        printf("Unsupported format: %s\n", ph->magic_nr);
        free(ph->magic_nr);
        return;
    }

	printf("Loaded %s\n", reading_file);
}

int vercoord(char *coordinates) {
    int8_t *p;
    int32_t how_many = 0;
	
    	p = (signed char*)strtok(coordinates, " ");
     while (p) {
        int e_number = true; 
        for (int i = 0; p[i] != '\0' && p[i] !='\n'; i++) {
            if (islower(p[i])) { 
                e_number = false;
                break;
            }
        }
        if (e_number == false) {
            break;
        }
        how_many++;
        p = (signed char*)strtok(NULL, " ");
    }

    return how_many == 4 ? true : false; 
}

void swap(int *a, int *b) {
    *a += *b;
    *b = *a - *b;
    *a -= *b;
}

void verlens(image **ph) {
    if ((*ph)->coord[0] > (*ph)->coord[2]) {
        swap((int*)&(*ph)->coord[0], (int*)&(*ph)->coord[2]);
    }

    if ((*ph)->coord[1] > (*ph)->coord[3]) {
        swap((int*)&(*ph)->coord[1], (int*)&(*ph)->coord[3]);
    }
}




void versel(image *ph, short *e_ok) {
    
    for (int i = 0; i < 4; i++) {
        if ((i % 2 == false && (ph->coord[i] > ph->coloane || (int)ph->coord[i] < 0)) || 
            (i % 2 == true && (ph->coord[i] > ph->rand || (int)ph->coord[i] < 0))) {
            *e_ok = false;
            return; 
        }
    }

int start_x = *(ph->coord);
int start_y = *(ph->coord + 1);
int end_x = *(ph->coord + 2);
int end_y = *(ph->coord + 3);

int is_same_point = (start_x == end_x && start_y == end_y);
int is_out_of_bounds = (start_x >= (int)ph->coloane || start_y >= (int)ph->rand);

if (is_same_point || is_out_of_bounds) {
    *e_ok = 0;  
    return;
}
}


void select_pixels(image *ph, char *given_img) {
    if (ph->state.e_loaded == false) {
        printf("No image loaded\n");
        return;
    }

  
    char *copy_of_coord = strdup(given_img);
    int e_ok = vercoord(copy_of_coord);
    free(copy_of_coord);

    if (e_ok == false) {
        printf("Invalid command\n");
        return;
    }

   
    if (ph->first_selection == false) ph->first_selection = true;

    if (ph->selection_counter == false && ph->state.e_first == false) {
        ph->coord = malloc(4 * sizeof(int));
        ph->state.e_first = true;
    }

    
    int *old_coord = NULL;
    if (ph->selection_counter) {
        old_coord = (int *)malloc(4 * sizeof(int));
        memcpy(old_coord, ph->coord, 4 * sizeof(int));
    }

  
    int index = 0;
    char *tmp = strtok(given_img, " ");
    while (tmp && index < 4) {
        ph->coord[index++] = atoi(tmp);
        tmp = strtok(NULL, " ");
    }

   
    verlens(&ph);
    versel(ph, (short int*)&e_ok);

    if (e_ok == true) {
        printf("Selected %d %d %d %d\n", ph->coord[0], ph->coord[1], ph->coord[2], ph->coord[3]);
        
        int e_full_selection = 
            (ph->coord[2] - ph->coord[0] == ph->coloane) && 
            (ph->coord[3] - ph->coord[1] == ph->rand);

        ph->selection_counter = e_full_selection ? false : true;
        ph->state.e_all = e_full_selection ? true : false;
    } else {
        printf("Invalid set of coordinates\n");
        if (old_coord) {
            memcpy(ph->coord, old_coord, 4 * sizeof(int));
        }
    }

    free(old_coord);
}


void select_all(image *ph, int ok) {
    if (!ok) {
        ph->state.e_all = false;
        return;
    }

    if (!ph->state.e_loaded) {
        printf("No image loaded\n");
        return;
    }

    ph->state.e_all = true;
    ph->selection_counter = false;
    printf("Selected ALL\n");
}


void rotate(double ***mat, int x1, void *x2, int y1, void *y2, int e_all)
{
	int valx2 = *(int*)x2;
	int valy2 = *(int*)y2;
    int32_t nr_rand = valy2 - y1;
    int32_t nr_coloane = valx2 - x1;

    if (e_all == false && nr_rand != nr_coloane) {
        printf("The selection must be square\n");
        return;
    }

    int32_t min_row = e_all ? x1 : y1;
    int32_t max_row = e_all ? valx2 : valy2;
    int32_t min_col = e_all ? y1 : x1;
    int32_t max_col = e_all ? valy2 : valx2;

    int32_t total_elements = nr_rand * nr_coloane;
    double *temp_buffer = calloc(total_elements, sizeof(double));

    if (!temp_buffer) {
        printf("Memory allocation failed\n");
        return;
    }

    if (nr_rand < nr_coloane) {
        *mat = (double **)realloc(*mat, nr_coloane * sizeof(double *));
        for (int i = nr_rand; i < nr_coloane; i++) {
            (*mat)[i] = malloc(nr_coloane * sizeof(double));
        }
    } else if (e_all == true && nr_rand != nr_coloane) {
        for (int i = 0; i < nr_rand; i++) {
            (*mat)[i] = realloc((*mat)[i], nr_rand * sizeof(double));
        }
    }

    int32_t idx = 0;
    for (int i = valx2 - 1; i >= x1; i--) {
        for (int j = y1; j < valy2; j++) {
            memcpy(&temp_buffer[idx], &((*mat)[j][i]), sizeof((*mat)[j][i]));
			idx += 1;
        }
    }

    idx = 0;
    for (int i = min_row; i < max_row; i++) {
        for (int j = min_col; j < max_col; j++) {
            memcpy(&(*mat)[i][j], &temp_buffer[idx], sizeof((*mat)[i][j]));
			idx++;

        }
    }

    free(temp_buffer);

    if (e_all == true) {
        if (nr_rand > nr_coloane) {
            for (int i = nr_coloane; i < nr_rand; i++) {
                free((*mat)[i]);
            }
            *mat = realloc(*mat, nr_coloane * sizeof(double *));
        } else if (nr_rand < nr_coloane) {
            for (int i = 0; i < nr_coloane; i++) {
                (*mat)[i] = realloc((*mat)[i], nr_rand * sizeof(double));
            }
        }
    }
}


void rotate_all_colors(image *ph, int start_w, void *width, int start_h, void *height, int direction) {
    rotate(&ph->mat_r, start_w, width, start_h, height, direction);
    rotate(&ph->mat_g, start_w, width, start_h, height, direction);
    rotate(&ph->mat_b, start_w, width, start_h, height, direction);
}

void rotimg(image *ph, char *string_angle)
{
    int angle = atoi(string_angle);

    if (!ph->state.e_loaded) {
        printf("No image loaded\n");
        return;
    }

    if (angle % 90 != 0 || angle > 360 || angle < -360) {
        printf("Unsupported rotation angle\n");
        return;
    }

    int effective_angle = (angle < 0) ? -angle : (360 - angle) % 360;

    if (!(ph->selection_counter || ph->state.e_all)) {
        printf("Rotated %d\n", angle);
        return;
    }

    while (effective_angle > 0) {
        if (ph->state.e_all) {
            uint32_t width = ph->coloane;
            uint32_t height = ph->rand;

            if (ph->state.e_rgb) {
                rotate_all_colors(ph, 0, &width, 0, &height, 1);
            } else {
                (void)rotate(&ph->img, 0, &width, 0, &height, 1);
            }

            ph->coloane = height;
            ph->rand = width;
        } else {
            int coords[4] = {ph->coord[0], ph->coord[2], ph->coord[1], ph->coord[3]};

            if (ph->state.e_rgb) {
                for (int i = 0; i < 3; ++i) {
                    double ***mat = (i == 0) ? &ph->mat_r : (i == 1) ? &ph->mat_g : &ph->mat_b;
                    rotate(mat, coords[0], &coords[1], coords[2], &coords[3], 0);
                }
            } else {
                rotate(&ph->img, coords[0], &coords[1], coords[2], &coords[3], 0);
            }
        }

        effective_angle -= 90;
    }
	fflush(stdout);
    fprintf(stdout,"Rotated %d\n", angle);
}

void eqcoords(int *x1, int *x2, int *y1, int *y2, image *ph)
{
    *x1 = ph->state.e_all ? 0 : ph->coord[0];
    *y1 = ph->state.e_all ? 0 : ph->coord[1];
    *x2 = ph->state.e_all ? ph->coloane : ph->coord[2];
    *y2 = ph->state.e_all ? ph->rand : ph->coord[3];
}

void crop(image *ph)
{
    if (!ph->state.e_loaded) {
        printf("No image loaded\n");
        return;
    }

    if (ph->selection_counter == false)
        return;

    int x1 = 0, x2 = 0, y1 = 0, y2 = 0;
    eqcoords(&x1, &x2, &y1, &y2, ph);

    int nr_coloane = 0;
    int nr_rand = 0;
    double *tmp_r = NULL, *tmp_g = NULL, *tmp_b = NULL, *tmp = NULL;
	nr_coloane += x2 - x1;
	nr_rand += y2 - y1;
    if (ph->state.e_rgb) {
        tmp_r = malloc(nr_rand * nr_coloane * sizeof(double));
        tmp_g = malloc(nr_rand * nr_coloane * sizeof(double));
        tmp_b = malloc(nr_rand * nr_coloane * sizeof(double));
    } else {
        tmp = malloc(nr_rand * nr_coloane * sizeof(double));
    }
	
    int index = false;
    for (int i = y1; i < y2; i++) {
        for (int j = x1; j < x2; j++) {
            if (ph->state.e_rgb) {
            memcpy(&tmp_r[index], &ph->mat_r[i][j], sizeof(ph->mat_r[i][j]));
			memcpy(&tmp_g[index], &ph->mat_g[i][j], sizeof(ph->mat_g[i][j]));
			memcpy(&tmp_b[index], &ph->mat_b[i][j], sizeof(ph->mat_b[i][j]));
			index++;

            } else {
            memcpy(&tmp[index], &ph->img[i][j], sizeof(ph->img[i][j]));
			index++;

            }
        }
    }
	
    for (int i = nr_rand; i < (int)ph->rand; i++) {
        if (ph->state.e_rgb == false) {
           free(ph->img[i]);
        } else {
            free(ph->mat_r[i]);
            free(ph->mat_g[i]);
            free(ph->mat_b[i]);
        }
    }

    if (ph->state.e_rgb) {
        ph->mat_r = realloc(ph->mat_r, nr_rand * sizeof(double *));
        ph->mat_g = realloc(ph->mat_g, nr_rand * sizeof(double *));
        ph->mat_b = realloc(ph->mat_b, nr_rand * sizeof(double *));
    } else {
        ph->img = realloc(ph->img, nr_rand * sizeof(double *));
    }

    index = 0;
    for (int i = 0; i < nr_rand; i++) {
        for (int j = 0; j < nr_coloane; j++) {
            if (ph->state.e_rgb) {
			memcpy(&ph->mat_r[i][j], &tmp_r[index], sizeof(tmp_r[index]));
			memcpy(&ph->mat_g[i][j], &tmp_g[index], sizeof(tmp_g[index]));
			memcpy(&ph->mat_b[i][j], &tmp_b[index], sizeof(tmp_b[index]));
			index++;

            } else {
			memcpy(&ph->img[i][j], &tmp[index], sizeof(tmp[index]));
			index++;
            }
        }
    }
	
    if (ph->state.e_rgb == true) {
        free(tmp_r);
        free(tmp_g);
        free(tmp_b);
    } else {
        free(tmp);
    }

    ph->rand = nr_rand;
    ph->coloane = nr_coloane;
    ph->selection_counter = 0;
    ph->state.e_all = 1;

    printf("Image cropped\n");
}

int rotunjire(double x)
{
	return (int)(x + 0.5);
}

double triaza(double x)
{
    return (x < 0) ? 0 : (x > 255) ? 255 : x;
}

double filt_sum(double **mat, const char *filt_name, int row, int col, nucleu filt)
{	row++;
	col++;
    float sum = 0;
    int index = 0;

    for (int i = row - 2; i <= row; i++) {
        for (int j = col - 2; j <= col; j++,index++) {
            double weight = 0;

            if (strcmp(filt_name, "EDGE") == 0)
                weight = filt.edge[index];
            else if (strcmp(filt_name, "SHARPEN") == 0)
                weight = filt.sharpen[index];
            else if (strcmp(filt_name, "BLUR") == 0)
                weight = filt.blur[index] / 9.0;
            else if (strcmp(filt_name, "GAUSSIAN_BLUR") == 0)
                weight = filt.gauss[index] / 16.0;

            sum += mat[i][j] * weight;
        }
    }

    return triaza(sum);
}


int verify(const char *parameter)
{
    const char *filt_names = "EDGE SHARPEN BLUR GAUSSIAN_BLUR";
    return strstr(filt_names, parameter) != NULL;
}

void apply_filt(image *ph, const char *filt_name, nucleu filt)
{
    if (!ph->state.e_loaded) {
        printf("No image loaded\n");
        return;
    }

    if (!verify(filt_name)) {
        printf("APPLY parameter invalid\n");
        return;
    }

    if (strcmp(ph->magic_nr, "P2") == 0 || strcmp(ph->magic_nr, "P5") == 0) {
        printf("Easy, Charlie Chaplin\n");
        return;
    }

    int x1 =0 , y1 = 0, x2 = 0, y2 = 0;
    eqcoords(&x1, &x2, &y1, &y2, ph);
	int nr_coloane = 0;
    int nr_rand = 0;
    
    double *tmp_r = NULL, *tmp_g = NULL, *tmp_b = NULL, *tmp = NULL;
	nr_rand += y2 - y1;
    nr_coloane += x2 - x1;

    if (ph->state.e_rgb) {
        tmp_r = calloc(nr_rand * nr_coloane, sizeof(double));
        tmp_g = calloc(nr_rand * nr_coloane, sizeof(double));
        tmp_b = calloc(nr_rand * nr_coloane, sizeof(double));
    } else {
        tmp = calloc(nr_rand * nr_coloane, sizeof(double));
    }

    int index = false;
    for (int i = y1; i < nr_rand + y1; i++) {
        for (int j = x1; j < nr_coloane + x1; j++) {
            if (ph->state.e_rgb == true) {
                if (i > 0 && i < (int)ph->rand - 1 && j > 0 && j < (int)ph->coloane - 1) {

				double R = filt_sum(ph->mat_r, filt_name, i, j, filt);
				double G = filt_sum(ph->mat_g, filt_name, i, j, filt);
				double B = filt_sum(ph->mat_b, filt_name, i, j, filt);
                memcpy(&tmp_r[index], &R, sizeof(double));
				memcpy(&tmp_g[index], &G, sizeof(double));
				memcpy(&tmp_b[index], &B, sizeof(double));;
                } else {
                    memcpy(&tmp_r[index], &ph->mat_r[i][j], sizeof(tmp_r[index]));
					memcpy(&tmp_g[index], &ph->mat_g[i][j], sizeof(tmp_g[index]));
					memcpy(&tmp_b[index], &ph->mat_b[i][j], sizeof(tmp_b[index]));

                }
                index++;
            } else {
                tmp[index] = (i > 0 && i < (int)ph->rand - 1 && j > 0 && j < (int)ph->coloane - 1) ? filt_sum(ph->img, filt_name, i, j, filt) : ph->img[i][j];
                index++;
            }
        }
    }

    index = 0;
    for (int i = y1; i < y2; i++) {
        for (int j = x1; j < x2; j++) {
            if (ph->state.e_rgb) {
			memcpy(&ph->mat_r[i][j], &tmp_r[index], sizeof(tmp_r[index]));
			memcpy(&ph->mat_g[i][j], &tmp_g[index], sizeof(tmp_g[index]));
			memcpy(&ph->mat_b[i][j], &tmp_b[index], sizeof(tmp_b[index]));
			index++;

            } else {
			memcpy(&ph->img[i][j], &tmp[index], sizeof(tmp[index]));
			index++;
            }
        }
    }

    if (ph->state.e_rgb) {
        free(tmp_r);
        free(tmp_g);
        free(tmp_b);
    } else {
        free(tmp);
    }

    printf("APPLY %s done\n", filt_name);
}


void save(image *ph, const char *saving_file) {
    if (!ph->state.e_loaded) {
        printf("No image loaded\n");
        return;
    }

    FILE *f;
    int e_binary = strstr(saving_file, "ascii") == NULL;

    if (e_binary) {
        f = fopen(saving_file, "wb");
    } else {
        char *file_name = strtok((char *)saving_file, " ");
        f = fopen(file_name, "wt");
        saving_file = file_name;
    }

    if (!f) {
        printf("Error saving file\n");
        return;
    }

    printf("Saved %s\n", saving_file);

    if (!e_binary) {
        if (strcmp(ph->magic_nr, "P5") == 0) strcpy(ph->magic_nr, "P2");
        else if (strcmp(ph->magic_nr, "P6") == 0) strcpy(ph->magic_nr, "P3");
    } else {
        if (strcmp(ph->magic_nr, "P2") == 0) strcpy(ph->magic_nr, "P5");
        else if (strcmp(ph->magic_nr, "P3") == 0) strcpy(ph->magic_nr, "P6");
    }

    fprintf(f, "%s\n%d %d\n%d\n", ph->magic_nr, ph->coloane, ph->rand, ph->intensity);

    for (int i = 0; i < (int)ph->rand; i++) {
        for (int j = 0; j < (int)ph->coloane; j++) {
            if (ph->state.e_rgb) {
                unsigned char r = (unsigned char)rotunjire(ph->mat_r[i][j]);
                unsigned char g = (unsigned char)rotunjire(ph->mat_g[i][j]);
                unsigned char b = (unsigned char)rotunjire(ph->mat_b[i][j]);

                if (e_binary) {
                    fwrite(&r, sizeof( char), 1, f);
                    fwrite(&g, sizeof( char), 1, f);
                    fwrite(&b, sizeof( char), 1, f);
                } else {
                    fprintf(f, "%d %d %d ", r, g, b);
                }
            } else {
                unsigned char pixel = (unsigned char)rotunjire(ph->img[i][j]);

                if (e_binary) {
                    fwrite(&pixel, sizeof(unsigned char), 1, f);
                } else {
                    fprintf(f, "%d ", pixel);
                }
            }
        }

        if (!e_binary) {
            fprintf(f, "\n");
        }
    }

    fclose(f);
}

void close(image *ph, nucleu *filt)
{   if (ph->state.e_loaded) {
		if (ph->state.e_rgb) {
			FREE_2D_ARRAY(ph->mat_r, ph->rand);
			FREE_2D_ARRAY(ph->mat_g, ph->rand);
			FREE_2D_ARRAY(ph->mat_b, ph->rand);
		} else {
			FREE_2D_ARRAY(ph->img, ph->rand);
		}
		free(ph->magic_nr);
		if (ph->first_selection == true)
			free(ph->coord);
	}
	
void* filter_pointers[] = {filt->edge, filt->blur, filt->gauss, filt->sharpen};

free(filter_pointers[0]);
free(filter_pointers[1]);
free(filter_pointers[2]);
free(filter_pointers[3]);
}


void ack(int *len, char **cmd, char **operation, char **copy_of_cmd) {
    while (*len > 0 && isspace((*cmd)[*len])) {
        (*cmd)[*len] = '\0';
        (*len)--;
    }

    strcpy(*copy_of_cmd, *cmd);

    if ((*operation = strtok(*cmd, " "))) {
        *len = strlen(*operation);
    }
}

char *alloc_char(char **cmd, int len, char *copy_of_cmd)
{	memcpy(*cmd, copy_of_cmd, strlen(copy_of_cmd) + 1);


	char *tmp = malloc(BASIC_SIZE);
	strcpy(tmp, *cmd + len + 1);

	return tmp;
}


void handle_load(image *ph, char **cmd, int len, char *copy_of_cmd) {
    char *tmp = alloc_char(cmd, len, copy_of_cmd);
    select_all(ph, 0);
    incaraca(ph, tmp);
    free(tmp);
}

void handle_select(image *ph, char **cmd, int len, char *copy_of_cmd) {
    char *tmp = alloc_char(cmd, len, copy_of_cmd);
    if (strcmp(copy_of_cmd, "SELECT ALL") == 0) {
        select_all(ph, 1);
    } else {
        select_pixels(ph, tmp);
    }
    free(tmp);
}

void handle_rotate(image *ph, char **cmd, int len, char *copy_of_cmd) {
    char *tmp = alloc_char(cmd, len, copy_of_cmd);
    rotimg(ph, tmp);
    free(tmp);
}

void handle_save(image *ph, char **cmd, int len, char *copy_of_cmd) {
    char *tmp = alloc_char(cmd, len, copy_of_cmd);
    save(ph, tmp);
    free(tmp);
}

void handle_apply(image *ph, char **cmd, int len, char *copy_of_cmd, nucleu filt) {
    char *tmp = alloc_char(cmd, len, copy_of_cmd);
    int apply_len = strlen(*cmd);

    if (ph->state.e_loaded) {
        if (apply_len > 5) {
            apply_filt(ph, tmp, filt);
        } else {
            printf("Invalid command\n");
        }
    } else {
        printf("No image loaded\n");
    }
    free(tmp);
}

void handle_exit(image *ph, nucleu *filt, char *cmd, char *copy_of_cmd) {
    close(ph, filt);
    free(cmd);
    free(copy_of_cmd);

    if (ph->state.e_loaded == 0) {
        printf("No image loaded\n");
    }
}

void process_command(image *ph, nucleu *filt, char **cmd, char *copy_of_cmd, char *operation, int len) {
    if (strcmp(operation, "LOAD") == 0) {
        handle_load(ph, cmd, len, copy_of_cmd);
    } else if (strcmp(operation, "SELECT") == 0) {
        handle_select(ph, cmd, len, copy_of_cmd);
    } else if (strcmp(operation, "ROTATE") == 0) {
        handle_rotate(ph, cmd, len, copy_of_cmd);
    } else if (strcmp(*cmd, "CROP") == 0) {
        if (ph->selection_counter) {
            crop(ph);
        } else {
            if (ph->state.e_loaded == 0) {
                printf("No image loaded\n");
            } else {
                printf("Image cropped\n");
            }
        }
    } else if (strcmp(operation, "SAVE") == 0) {
        handle_save(ph, cmd, len, copy_of_cmd);
    } else if (strcmp(operation, "EXIT") == 0) {
        handle_exit(ph, filt, *cmd, copy_of_cmd);
        exit(0);
    } else if (strcmp(operation, "APPLY") == 0) {
        handle_apply(ph, cmd, len, copy_of_cmd, *filt);
    } else if (strcmp(operation, "EQUALIZE") == 0) {
        if (ph->state.e_loaded == false) {
            printf("No image loaded\n");
        } else if (ph->state.e_rgb == true) {
            printf("Black and white image needed\n");
        } else {
            printf("Equalize done\n");
        }
    } else if (strcmp(operation, "HISTOGRAM") == 0) {
        if (ph->state.e_loaded == false) {
            printf("No image loaded\n");
        } else if (ph->state.e_rgb == true) {
            printf("Black and white image needed\n");
        } else if ((strcmp(ph->magic_nr, "P2") == 0) || strcmp(ph->magic_nr, "P5") == 0) {
            printf("Invalid command\n");
        } else {
            printf("Equalize done\n");
        }
    } else {
        printf("Invalid command\n");
    }
}

int main(void) {
    nucleu filt;
    image ph;

    char *cmd = (char *)malloc(BASIC_SIZE);
    char *copy_of_cmd = (char *)malloc(BASIC_SIZE);

    initializare_filt(&filt);
    initializare(&ph);
    do {
        char *operation;
        fgets(cmd, BASIC_SIZE, stdin);

        if (strcmp(cmd, "\n") == false) {
            continue;
        }

        int len = strlen(cmd) - 1;
        ack(&len, &cmd, &operation, &copy_of_cmd);

        process_command(&ph, &filt, &cmd, copy_of_cmd, operation, len);
    } while(1);

    return 0;
}
