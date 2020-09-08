# include <stdlib.h> //standard inclusion
# include <stdio.h> //standard inclusion
# include <math.h> //for cos/sin operations
//# include <time.h>
//# include <omp.h>
# include <dirent.h> //for directory surfing
//# include <sys/stat.h> 
//# include <sys/types.h>
# include "julia_gen.h" /*custom header file for Julia Set generator...
reference is given in julia_gen.c as it is adapted from an existing...
mandelbrot openmp code*/

//gcc julia_gen.c -c -fopenmp -lm
//gcc julia_gen.o Mid_mitcha14.c -o Mid_mitcha14 -fopenmp -lm


int mkdir(const char *pathname, mode_t mode); //prototype for new folder

typedef struct { 
	unsigned char r,g,b; //RGB values (0:255) corresponding to pixel in image
} color; //datatype for pixels comprising any '.ppm' formated image

typedef struct {
	int m, n; //dimensions of image (pixels x pixels)
	int count; /*determines scaling of color map for iterated functions...
	maximum iteration divergence threshold before being ubiquitously 255 */
	double rotation; //angle in Euler's formula for Julia complex parameter
	double scale; //amplitute in Euler's formula for Julia complex parameter
	char *file_ext; /*...
	string with file extension of .ppm file created from call to img_gen()*/
	color *data; //unrolled array of pixels each with color values RGB
} img; //Julia Set image datatype linked to outputs of julia_gen.c

typedef struct {
	double *min; //minimum of varied parameter in train set (scale/rotation)
	double *max; //maximum of varied parameter in train set (scale/rotation)
	int num_dir; //total number of training sets (directories)
	int *num_image; //array of image quantities in each training set
	int *scale_flags;
	char **dir_ext; 
	/*array of strings with the extensions to directories corresponding to... 
		training sets generated with set_gen during present or previous runs*/
	img **images; /*2d array containing training sets in first dimension and...
	set images in second dimension.  Index of dimension one corresponds to...
	that of num_image and dir_ext*/
} dir; //Current directory, containing any number of training sets

static int num_samples; //number of samples compared to each training set


char * set_gen(int m, int n,         int count_max,    int num_image, 
					  double scale,  double rotation,  int scale_flag,
					  int midp_flag, double min_range, double max_range)
/* Purpose of function is to generate a training set of Julia set images...
by varing the complex parameter through a scaling or rotation (scale_flag)...
in Euler's formula between specified minimum and maximum bounds...
determined with either explicit declaration or a midpoint and perturbation...
(midp_flag)

Additional inputs include the image dimensions m, n...
a maximum iteration number count_max which decides when a point is...
bounded and hence an interior point of the Julia set with color r=g=b=255

Outputs a string with the extension to the training set directory...
just created, not actually used because of the parsing function grab_sets()*/
{
  double a; // variable representing the axis (scale/rotation) to increment
  double c1, c2; // real and imaginary parts of complext parameter
  char *output_dir = malloc(64); //training set directory
  char *output_filename = malloc(128); //training set image directory
  if (scale_flag == 1) //i.e. if scale is being incremented
    sprintf(output_dir, "./%d of %dmax %dx%d, a*exp(%lf*i) [%lf,%lf]/", 
  			num_image, count_max, m, n, rotation,min_range, max_range);
  else //i.e. if rotation is being incremented
    sprintf(output_dir, "./%d of %dmax %dx%d, %lf*exp(a*i) [%lf,%lf]/", 
  			num_image, count_max, m, n, scale, min_range, max_range); /*...
  Writes the name of the directory to contain information about the set...
  as a very specific format meant to be read and interpreted disjointly...
  from the following image generation. Custom image/set naming not supported*/
  mkdir(output_dir, 0777); //new folder with formatted name
  for (int in = 0; in < num_image; in++) { //iterating over image count in set
    a = min_range + in * (max_range - min_range) / (num_image-1);/*...
    setting value of variable incremented between max_range and min_range...
    each value corresponds to its own image*/
    sprintf(output_filename, "%s%f.ppm", output_dir, a); /*...
    set image extension to be the value of 'a' and file type to be '.ppm'*/

    if (scale_flag == 0) {// i.e. if rotation is being incremented
      c1 = scale*cos(a); //real part of Euler's formula, rotation variable
      c2 = scale*sin(a); //imag part of Euler's formula, rotation variable
    } else { // i.e. if scale is being incremeneted
      c1 = a*cos(rotation); //real part of Euler's formula, scale variable
      c2 = a*sin(rotation); //imag part of Euler's formula, scale variable
    }
    img_gen(m, n, count_max, scale, rotation, c1, c2, output_filename); /*
    call to function in 'julia_gen.h' that produces a Julia set with the...
    provied parameters*/
  }
  free(output_filename); //free malloc
  return output_dir; //again, will not be useful since grab_sets() obtains all
}


static dir *grab_sets(const char *input_dir) /*...
Purpose of function is to parse the current directory (input_dir) for info on...
sub-directories with the proper naming convention indicating they contain...
a unique set of training images produced by set_gen()

Parameters that determined the generation of the sets and containing images..
as well as the appropiate file extensions are stored for later use...
note that the pixel data from the image is not read in at this time

Output is a struct of datatype dir for the current directory*/
{
	dir *directories; //structure of current directory, final output
	img **images; //2d array for contained training sets and associated images
    struct dirent *entry, *sub_entry; //from <dirent.h> for entry in directory
    DIR *working = opendir(input_dir); //from <dirent.h>, opens directory
    if (working == NULL){ //opendir returns NULL if couldn't open directory 
        printf("Could not open current directory" ); //error handling
        exit(1); 
    } 
    /* https://www.geeksforgeeks.org/c-program-list-files-sub-directories-directory/
	for reference used to guided listing of all files and sub-directories */

	int i, j; //index of training sets i, and index of contained images j
    int *num_image = (int *)malloc(200 * sizeof(int)); /*preallocate array...
    of image counts in each directory with number of directories no more...
    than 200.  Extra space allocated since number of directories not yet...
    known and would otherwise require an additional loop through files*/
    int temp_num; //temporary used only to accept values in sscanf
    int num_dir = 0; //number of directories
    int count; //maximum iteration threshold
    int m, n; //dimensions
    double temp_min, temp_max; //temporary used only to accept values in sscanf
    int scale_flag, empty_flag; /* See previous for scale_flag, empty_flag is...
    active to skip image of training set not matching naming convention*/
    double scale, rotation; //See previous for scale/rotation
	double min, max; //See previous for min/max

    while ((entry = readdir(working)) != NULL) { // loop entries of current
  		if ((sscanf(entry->d_name, "%d of %dmax %dx%d, %lf*exp(a*i) [%lf,%lf]", 
  					&temp_num, &count, &m, &n, &scale, 
  					&temp_min, &temp_max) == 7)
  		 || (sscanf(entry->d_name, "%d of %dmax %dx%d, a*exp(%lf*i) [%lf,%lf]", 
					&temp_num, &count, &m, &n, &rotation, 
					&temp_min, &temp_max) == 7)) { //if entry matches naming

  			num_image[num_dir] = 0; //initialize image count of directory
  			DIR *sub_working = opendir(entry->d_name);//open directory

	  		while ((sub_entry = readdir(sub_working)) != NULL) {//loop images
	  			if (sscanf(sub_entry->d_name, "%lf.ppm", &rotation) == 1) {
	  			//if image matches naming
	  				num_image[num_dir]++; //increment up image count
	  			}
	  		}
	  		num_dir++; //increment up directory count
	  	}
  	};
  	closedir(working);//close directory to reset while loop

  	working = opendir(input_dir); //reopen at start
  	directories = (dir *)malloc(sizeof(dir));//preallocation for struct
  	directories->num_dir = num_dir; //number of training sets in current
  	directories->images = (img **)malloc(num_dir * sizeof(img *));
  	directories->min = (double *)malloc(num_dir * sizeof(double));
  	directories->max = (double *)malloc(num_dir * sizeof(double));
  	directories->num_image = (int *)malloc(num_dir * sizeof(int));
  	directories->scale_flags = (int *)malloc(num_dir * sizeof(int));
  	directories->dir_ext = malloc (num_dir * sizeof(char *));
  	/*Preallocation for an array of directories of size num_dir through...
  	corresponding arrays of various struct dir elements*/

  	i = 0;//Initialize directory increment
    while ((entry = readdir(working)) != NULL) { // loop entries of current
  		if (sscanf(entry->d_name, "%d of %dmax %dx%d, %lf*exp(a*i) [%lf,%lf]", 
  				   &temp_num, &count, &m, &n, &scale, 
  				   &temp_min, &temp_max) == 7) //if entry matches naming
  			scale_flag = 0; //indicates rotation has been incremented
  		else if (sscanf(entry->d_name, "%d of %dmax %dx%d, a*exp(%lf*i) [%lf,%lf]", 
  									   &temp_num, &count, &m, &n, &rotation, 
  									   &temp_min, &temp_max) == 7)
  			scale_flag = 1; //indicates scale has been incremented
  		else
  			scale_flag = 2; //indicates naming convention was not matched
  		if (scale_flag != 2) { //i.e. if entry is a valid training set
  			j = 0; //initialize image increment
  			directories->dir_ext[i] = malloc(64);//preallocate string for ext
			sprintf(directories->dir_ext[i], "./%s/", entry->d_name);//set ext
	  		directories->images[i] = (img *)malloc(num_image[i] * sizeof(img));
	  		//Preallocate num_image[i] arrays of images for each training set i
	  		directories->num_image[i] = num_image[i]; //retain for later use
	  		directories->scale_flags[i] = scale_flag; //retain for later use
	  		DIR *sub_working = opendir(entry->d_name); //open each training set
	  		while ((sub_entry = readdir(sub_working)) != NULL) {//loop images
	 			empty_flag = 0;//start empty, activate only if confirmed image
	  			if (scale_flag == 0) { // i.e. if rotation is incremeneted
	  				if (sscanf(sub_entry->d_name, "%lf.ppm", &rotation) == 1) 
	  				{ //if naming convention is matched for image
	  					if (j == 0) { //initializing minimum
	  						min = rotation;
	  						max = rotation;
	  					} else if (rotation < min) { //check if jth is lower
	  						min = rotation;
	  					} else if (rotation > max) { //check if jth is higher
	  						max = rotation;
	  					}
						empty_flag = 1; //confirmed image to proceed
					}
				} else if (sscanf(sub_entry->d_name, "%lf.ppm", &scale) == 1) {
	  					if (j == 0) { //initializing minimum
	  						min = scale;
	  						max = scale;
	  					} else if (scale < min) { //check if jth is lower
	  						min = scale;
	  					} else if (scale > max) { //check if jth is higher
	  						max = scale;
	  					}
						empty_flag = 1; //confirmed image to proceed
				}
				if (empty_flag == 1) { //proceeding with confirmed image
	  				directories->images[i][j].file_ext = malloc(128);
	  				//preallocating string for image extension
	  				sprintf(directories->images[i][j].file_ext, "./%s/%s", 
	  						entry->d_name, sub_entry->d_name);
	  				//write image extension using directory and image name

	  				directories->images[i][j].count = count;
	  				directories->images[i][j].m = m;
	  				directories->images[i][j].n = n;
	  				directories->images[i][j].rotation = rotation;
	  				directories->images[i][j].scale = scale;
	  				/*Assignment of struct img elements corresponding to...
	  				training set i and image j*/
	  				j++; //increment images
		  		}
			}
			closedir(sub_working);//close training set i
			if ((min != temp_min) || (max != temp_max) 
			 	 				  || (num_image[i] != temp_num)) {
      			fprintf(stderr, "Incorrect image number (missing or extra).\n");
			} /*non-critical error reporting for if range or number of...
			images found in a training set does not match that indicated...
			in the naming of the directory*/
			directories->min[i] = min;
			directories->max[i] = max;
			//final assignment of directory elements
			i++; //increment directory
		}		
	}
    closedir(working); //close current directory
  	free(num_image); //free malloc
    return directories; //returns dir struct directories missing data element
}

static img *train_sets(dir *directories) /*...
Purpose of function is to average color values of corresponding coordinates...
in the images of each training set.  Here, input and output for '.ppm'...
files is handled to read the images in input dir directories populated...
using grab_sets() and containing all information but the color data.

For each training set, images are read in one at a time to compute and...
save as a Portable PixelMap. Output is an array of averaged, trained...
images of datatype struct img to then compare samples against.

https://stackoverflow.com/questions/2693631/read-ppm-file-and-store-it-in-an-array-coded-with-c
for reference on reading PPM file and storing it as an array...
note here it is adapted for P3 format instead of P6*/
{
	img *trained_images; //array of trained images to be output
	trained_images = (img *)malloc(directories->num_dir * sizeof(img));
	//preallocation of trained images array of size num_dir
	double *temp_r;//temporary array of r values used for average
	double *temp_g;//temporary array of g values used for average
	double *temp_b;//temporary array of b values used for average
	char temp_buff[16];
	int temp_color;//temporary color used for error handling of image read
	int temp_m;//temporary dimension used for error handling of image read
	int temp_n;//temporary dimension used for error handling of image read
	mkdir("./Trained_images/", 0777);//create new folder trained images
	for (int i = 0; i < directories->num_dir; i++) {//iterate through dirs
		for (int j = 0; j < directories->num_image[i]; j++) {//iterate images
	  		FILE *fp; //to store file pointer as per <stdio.h>
			fp = fopen(directories->images[i][j].file_ext, "rb");//open image j
			
			/*Below is error handling to open .ppm image, see reference...
			note that these checks will not fail under normal operation...
			given that sets and images are guaranteed to be properly...
			formated and openable upon the successful completion of set_gen()...
			Each of below additionally serves to progress the scan beyond...
			the header and onto the array of RGB values.

			Same handling is duplicated further along and would ideally be...
			done through multiple calls to single function that would fully...
			handle file opening/closing and add data to dir directories*/
			if (!fp) {
      			fprintf(stderr, "Unable to open file '%s'\n", 
      							directories->images[i][j].file_ext);
      			exit(1);
 			} 
 			if (!fgets(temp_buff, sizeof(temp_buff), fp)) {
 				perror(directories->images[i][j].file_ext);
 				exit(1);
 			}
 			if (temp_buff[0] != 'P' || temp_buff[1] != '3') {
 				fprintf(stderr, "Invalid image format (must be 'P3')\n");
 				exit(1);
 			}
 			if (fscanf(fp, "%d  %d", &temp_m, &temp_n) != 2) {
 				fprintf(stderr, "Invalid image size (error loading '%s')\n",
 								directories->images[i][j].file_ext);
 			}
 			if (fscanf(fp, "%d", &temp_color) != 1) {
 				fprintf(stderr, "Invalid rgb component (error loading '%s')\n",
 								directories->images[i][j].file_ext);
 				exit(1);
 			}
 			while (fgetc(fp) != '\n'); //progress any additional blanks
 			
 			directories->images[i][j].data = (color *)malloc(
											 directories->images[i][j].m * 
											 directories->images[i][j].n * 
											 sizeof(color)); /*...
											 Preallocate unrolled array of...
											 pixels in image j, training set i*/
 			for (int k = 0; k < directories->images[i][j].m * 
 								directories->images[i][j].n; k++) {
 				//iterate through each k corresponding to unique image pixel
	  			fscanf(fp, "  %hhd  %hhd  %hhd", 
	  				   &directories->images[i][j].data[k].r,
	  				   &directories->images[i][j].data[k].g,
	  				   &directories->images[i][j].data[k].b); /*...
	  				   Progress three RGB values into the .ppm file*/

	  		}
	  		fclose(fp);//close image, having stored RGB for each pixel
	  		if (j == 0) { //initialize temp color values to be averaged
	  			temp_r = (double *)malloc(directories->images[i][j].m 
	  									  * directories->images[i][j].n
	  									  * sizeof(double));
				temp_g = (double *)malloc(directories->images[i][j].m
										  * directories->images[i][j].n
										  * sizeof(double));
				temp_b = (double *)malloc(directories->images[i][j].m
										  * directories->images[i][j].n
										  * sizeof(double));
				//preallocate unrolled arrays of pixel colors to be averaged

	  			trained_images[i].m = directories->images[i][j].m;
	  			trained_images[i].n = directories->images[i][j].n;
	  			trained_images[i].count = directories->images[i][j].count;
	  			//record parameters for future use

	  			for (int k = 0; k < trained_images[i].m
	  				 			    * trained_images[i].n; k++)  
  				{/*initialize averaged image pixel data with assignment from...
  				normalized pixel colors found in the j=0th image
				*/
	  				temp_r[k] = (double) directories->images[i][j].data[k].r
	  											/ directories->num_image[i];
	  				temp_g[k] = (double) directories->images[i][j].data[k].g
	  											/ directories->num_image[i];
	  				temp_b[k] = (double) directories->images[i][j].data[k].b
	  											/ directories->num_image[i];
	  			//normalized by total number of images in trained set
	  			}
	  		} else {//i.e. if no longer the first image in set
	  			for (int k = 0; k < trained_images[i].m
	  				 			    * trained_images[i].n; k++) { //loop pixels
  					temp_r[k] = temp_r[k]
  							    + (double)directories->images[i][j].data[k].r 
  								/ directories->num_image[i];
					temp_g[k] = temp_g[k]
								+ (double)directories->images[i][j].data[k].g
								/ directories->num_image[i];
  					temp_b[k] = temp_b[k]
  								+ (double)directories->images[i][j].data[k].b
  								/ directories->num_image[i];
  					/*summing cumulative averaged parts with kth images...
  					normalized pixel colors*/
	  			}
	  		} //done images in trained set i
	  		free(directories->images[i][j].data); //free large data array
  			free(directories->images[i][j].file_ext); //free image extension
	  	}
	  	trained_images[i].file_ext = malloc(128);
	  	//preallocate trained image file extension

	  	if (directories->scale_flags[i] == 0) { //i.e. if rotation incremented
	  		sprintf(trained_images[i].file_ext, 
	  				"./Trained_images/%d of %dmax %dx%d, %lf*exp(a*i) [%f,%f].ppm",
	  				directories->num_image[i], trained_images[i].count, 
	  				trained_images[i].m, trained_images[i].n, 
	  				directories->images[i][0].scale,
	  				directories->min[i],directories->max[i]);
	  	} else { //i.e. if scale incremented in training set
	  		sprintf(trained_images[i].file_ext, 
	  				"./Trained_images/%d of %dmax %dx%d, a*exp(%lf*i) [%f,%f].ppm",
	  				directories->num_image[i], trained_images[i].count, 
	  				trained_images[i].m, trained_images[i].n, 
	  				directories->images[i][0].rotation,
	  				directories->min[i],directories->max[i]);
	  	}/*Write file extension of trained image to contain same information...
	  	as the training set it was crafted from (long name is needed...
	  	to avoid writing over the same image with differing sets)*/ 

  		FILE *output_unit; //as per <stdio.h> for handling files
  		output_unit = fopen ( trained_images[i].file_ext, "wt" );
  		//create .ppm file with the appropiate naming convention
		fprintf ( output_unit, "P3\n" ); //PPM "magic number" formatting
		fprintf ( output_unit, "%d  %d\n", trained_images[i].m, 
										   trained_images[i].n );//dimensions
		fprintf ( output_unit, "%d\n", 255 ); //RGB component type
		for ( int k = 0; k < trained_images[i].m; k++ ) {
    		for ( int l = 0; l < trained_images[i].n; l++) {
        		fprintf ( output_unit, "  %d %d %d", 
        			(int) temp_r[k * trained_images[i].m + l],
        			(int) temp_g[k * trained_images[i].m + l],
        			(int) temp_b[k * trained_images[i].m + l]);
        		/*packing data array of averaged pixel colors in temporary...
        		arrays with unrolled indexes to .ppm file body*/
  			}
  			fprintf ( output_unit, "\n" );//useful for format, ignored in PPM
		}//iterating over rows and coloumns
		fclose(output_unit);//close trained image file
		printf("Averaged %s\n", trained_images[i].file_ext);
		//informing the user trained image was completed
	}
	return trained_images; /*Array of datatype struct img with size num_dir...
	Note that unlike training sets in grab_sets() and samples in...
	grab_samples(), the trained_images struct does not parse its directory...
	to find those created in previous runs. Only training sets present...
	during runtime are averaged and compared to samples. Note .data...
	elements are to be loaded in one at a time in the compare stage*/
}

static img *grab_samples(const char *input_dir)
{/*Purpose of function is to return array of datatype struct img containing...
all images present in the ./Samples/ directory.  Note .data elements are to...
be loaded in one at a time during the compare stage*/
    struct dirent *entry; //from <dirent.h> for entry in directory
	int i; //sample index
    int count; //maximum iteration threshold, see earlier
    int m, n; //sample image dimensions
    double scale, rotation; //See earlier for scale/rotation
	DIR *working = opendir(input_dir);//from <dirent.h>, opens directory
	img *images; //array for datatype struct img for sample images
	
	num_samples = 0; //intialize sample count
    while ((entry = readdir(working)) != NULL) {//loop through entries
  		if (sscanf(entry->d_name, "%dmax %dx%d, %lf*exp(%lf*i)", 
  						&count, &m, &n, &scale, &rotation) == 5)
	  		num_samples++;//increment up sample count upon matched naming
  	}
  	closedir(working);//close Samples directory to reset while

  	working = opendir(input_dir); //reopen to start at beginning of directory
  	images = (img *)malloc(num_samples * sizeof(img)); 
  	//preallocate array of sample image datatype struct img of size num_samples 
	i = 0; //initialize sample index
    while ((entry = readdir(working)) != NULL) {//loop through entries
  		if (sscanf(entry->d_name, "%dmax %dx%d, %lf*exp(%lf*i)", 
  						&count, &m, &n, &scale, &rotation) == 5) {
  					//sample naming convention matched and parameters stored
  					images[i].count = count;
	  				images[i].m = m;
	  				images[i].n = n;
	  				images[i].scale = scale;
  					images[i].rotation = rotation;
  					//assign image parameters used to generate sample
  					images[i].file_ext = malloc(128);//preallocate file ext
	  				sprintf(images[i].file_ext, "./Samples/%s", entry->d_name);
	  				//store file extension for later use in pulling data
	  				i++; //increment up sample index
		}	
  	}
  	return images; //return array of datatype struct img for sample images
}


double compare(img trained_image, img sample)
{/*Purpose of function is to calculate the correlation coefficent between...
an image from the set of trained images and an image from the set of samples...
Value is -1 for perfectly anti-correlation and 1 for perfectly correlation...

Inputs are two singular img datatypes from the arrays trained_images and...
samples and are nearly entirely populated apart from the .data element 

See https://www.geeksforgeeks.org/program-find-correlation-coefficient/
for reference on the correlation coefficient calculation*/

	double correlation[3];/*array for the three color components containing...
	correlation coeffecients between trained_image and sample*/

	char temp_buff[16];//temporary buff for error handling, see earlier
	int temp_color;//temporary color used for error handling of image read
	int temp_m;//temporary dimension used for error handling of image read
	int temp_n;//temporary dimension used for error handling of image read

	int k, l; //pixel index k and color component index l
	double temp_trnsam[6]; //enables iterating through l with for loop
	double sum_trn[3] = {0.0,0.0,0.0}; //initialization of trained image sum
	double sum_sam[3] = {0.0,0.0,0.0}; //initialization of sample image sum
	double sum_trnsam[3] = {0.0,0.0,0.0};
	//initialization of trained image * sample image sum
	double sum_trn2[3] = {0.0,0.0,0.0};
	//initialization of trained image squared sum
	double sum_sam2[3] = {0.0,0.0,0.0};
	//initialization of trained image squared sum

	FILE *fp;//file pointer, see earlier
	fp = fopen(trained_image.file_ext, "rb");//opening trained image

	if (!fp) {
		fprintf(stderr, "Unable to open file '%s'\n", 
						trained_image.file_ext);
		exit(1);
	} 
	if (!fgets(temp_buff, sizeof(temp_buff), fp)) {
		perror(trained_image.file_ext);
		exit(1);
	}
	if (temp_buff[0] != 'P' || temp_buff[1] != '3') {
		fprintf(stderr, "Invalid image format (must be 'P3')\n");
		exit(1);
	}
	if (fscanf(fp, "%d  %d", &temp_m, &temp_n) != 2) {
		fprintf(stderr, "Invalid image size (error loading '%s')\n",
						trained_image.file_ext);
	}
	if (fscanf(fp, "%d", &temp_color) != 1) {
		fprintf(stderr, "Invalid rgb component (error loading '%s')\n",
						trained_image.file_ext);
		exit(1);
	}
	while (fgetc(fp) != '\n');
	//error handling and progressing beyond header, see earlier

	trained_image.data = (color *)malloc( trained_image.m * trained_image.n * 
										  sizeof(color));
	//preallocating .data element of trained image with datatype img
	for (int k = 0; k < trained_image.m * trained_image.n; k++) {
		fscanf(fp, "  %hhd  %hhd  %hhd", &trained_image.data[k].r,
									     &trained_image.data[k].g,
									     &trained_image.data[k].b);
		//obtaining RGB colors from PPM file body and storing to img.data
	} //iterate through pixels
	fclose(fp);//close trained image now that it is stored

	fp = fopen(sample.file_ext, "rb");//opening sample image
	if (!fp) {
		fprintf(stderr, "Unable to open file '%s'\n", 
						sample.file_ext);
		exit(1);
	} 
	if (!fgets(temp_buff, sizeof(temp_buff), fp)) {
		perror(sample.file_ext);
		exit(1);
	}
	if (temp_buff[0] != 'P' || temp_buff[1] != '3') {
		fprintf(stderr, "Invalid image format (must be 'P3')\n");
		exit(1);
	}
	if (fscanf(fp, "%d  %d", &temp_m, &temp_n) != 2) {
		fprintf(stderr, "Invalid image size (error loading '%s')\n",
						sample.file_ext);
	}
	if (fscanf(fp, "%d", &temp_color) != 1) {
		fprintf(stderr, "Invalid rgb component (error loading '%s')\n",
						sample.file_ext);
		exit(1);
	}
	while (fgetc(fp) != '\n');
	//error handling and progressing beyond header, see earlier

	sample.data = (color *)malloc( sample.m * sample.n * sizeof(color));
	//preallocating .data element of sample image with datatype img
	for (k = 0; k < sample.m * sample.n; k++) {//iterate through pixels
		fscanf(fp, "  %hhd  %hhd  %hhd", &sample.data[k].r,
									     &sample.data[k].g,
									     &sample.data[k].b);
		//obtaining RGB colors from PPM file body and storing to img.data
	}
	fclose(fp);//close sample image now that it is stored

	for (k = 0; k < sample.m * sample.n; k++) { //iterate through pixels
		//sample/trained image dimensions equal because of preselection
		temp_trnsam[0] = (double)trained_image.data[k].r;
		temp_trnsam[1] = (double)trained_image.data[k].g;
		temp_trnsam[2] = (double)trained_image.data[k].b;
		temp_trnsam[3] = (double)sample.data[k].r;
		temp_trnsam[4] = (double)sample.data[k].g;
		temp_trnsam[5] = (double)sample.data[k].b;
		//temporary array of color values to be able to use l for loop
		for (l = 0; l < 3; l++) {//iterate over color components (RGB)
			sum_trn[l] = sum_trn[l] + temp_trnsam[l] / 255.0;
			//summation of trained image pixels (l=0,1,2)
			sum_sam[l] = sum_sam[l] + temp_trnsam[l+3] / 255.0;
			//summation of sample image pixels (l=3,4,5)
			sum_trnsam[l] = sum_trnsam[l] + temp_trnsam[l]
							* temp_trnsam[l+3] / (255.0 * 255.0);
			//summation of trained image * sample image pixels
		  	sum_trn2[l] = sum_trn2[l] + temp_trnsam[l]
						  * temp_trnsam[l] / (255.0 * 255.0);
			//summation of trained image pixels squared
		    sum_sam2[l] = sum_sam2[l] + temp_trnsam[l+3]
						  * temp_trnsam[l+3] / (255.0 * 255.0);
			//summation of sample image pixels squared
		}	

	}
	for (l = 0; l < 3; l++) { //iterate over color components
		correlation[l] = (double)(sample.m * sample.n * sum_trnsam[l]
							   	  - sum_trn[l] * sum_sam[l]) 
							   	  / sqrt((sample.m * sample.n * sum_trn2[l] 
							   	  - sum_trn[l] * sum_trn[l])
							   	  * (sample.m * sample.n * sum_sam2[l]
							   	  - sum_sam[l] * sum_sam[l]));	
		//correlation coefficient calculation, see reference 
	}
	return (correlation[0] + correlation[1] + correlation[2])/3.0;
	//averages the correlation coefficient of each color component (RGB)
}

int main() 
{
	dir *directories; //array of datatype struct dir for each training set
	img *trained_images; //array of datatype struct img of trained images
	img *samples; //array of datatype struct img of sample images to compare

  	int scale_flag = 1; /*FALSE to increment rotation in Euler's formula...
  	for complex parameter in generation of training set of Julia set images...
  	TRUE to increment scale*/
  	int midp_flag = 0; /*FALSE to specify bounds of training set parameter...
  	increment through explicit assignment of min_range/max_range...
  	TRUE to specify bounds through explicit assignment of a midpoint and...
  	a perturbation with min_range/max_range calculated accordingly*/ 

  	int i, j; //index of training sets/trained images i and samples j

	int m, n; //dimensions of images in training set generation
	int count_max; //maximum iteration threshold for Julia set
	int num_image; //number of images to be generated in set_gen()
	double c1, c2; //Real and imaginary parts of complex parameter in Julia set
  	double scale, rotation; //in Euler's formula for complex parameter
 	double perturb, midpoint; //used to find min/max when midp_flag = 1

    double min_range, max_range; //min and max of incremented training set

	double correlation; //correlation coefficient returned from compare()

    char *output_dir = malloc(32); /*extension returned from set_gen,...
    unused afterwards since grab_sets() finds extensions of all training sets...
    in the current directory*/

/*=++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=*/
 	if (scale_flag == 0) {//i.e. if rotation is to be incremented
	    	perturb = 0.075; 		midpoint = M_PI*0.555;
 			//perturb = M_PI;			midpoint = 0.0;
	    	//Change these values when using midp_flag == 1 to specify min/max
	    if (midp_flag == 0) {//i.e. not using midpoint and perturbation method
	    	min_range = 0.95*M_PI; 	max_range = M_PI;
	    	//Change these values when using midp_flag == 0 for increment range

	    }
  	} else {//i.e. if scale is to be incremented
  			perturb = 0.25; 		midpoint = 0.7885;
  			//Change these values when using midp_flag == 1 to specify min/max
	    if (midp_flag == 0) {
	    	min_range = 0.65; 		max_range = 0.85;
	    	//Change these values when using midp_flag == 0 for increment range
	    }
	}
/*=++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=*/

	if (midp_flag == 1) {//i.e. using midpoint and perturbation method
	    min_range = midpoint - perturb; //assignment of minimum increment
	    max_range = midpoint + perturb; //assignment of maximum increment
	}

/*=++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=*/
	m = 500; n = 500; count_max = 60; num_image = 5;
	scale = 0.7885; rotation = 0.56*M_PI;
	/*Change these values to specify parameters of training set to be...
	generated. Note that scale_flag will determine whether either scale...
	or rotation will instead be incremented between min and max*/
/*=++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=*/
	
	output_dir = set_gen(m, n, count_max, num_image, scale, rotation, 
						 scale_flag, midp_flag, min_range, max_range);
	//call to function that generates a training set

	directories = grab_sets(".");
	/*Call to function to parse current directory for training sets...
	able to handle training sets generated in previous runs provided...
	the naming convention has not been altered*/
	printf("\n");
	trained_images = train_sets(directories);
	/*Call to function to average each training set into an array of...
	trained images. Note that trained images from previous runs are...
	not handled by this code without their corresponding training set*/
	printf("\n");
	
	mkdir("./Samples/", 0777);//create a directory to store sample images
	char *output_filename = malloc(128);//preallocate space for image ext

/*=++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=*/
	count_max = 60;
	scale = 0.7885; rotation = 0.555*M_PI;
	/*Change these values to specify parameters of a sample image to be...
	generated. Note that dimension is left the same as that to generate...
	the training set to guarantee at least one training set can be compared*/
/*=++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=*/

	c1 = scale*cos(rotation); c2 = scale*sin(rotation);
	//Real and imaginary parts Euler's formula for complex parameter
	sprintf(output_filename, "./Samples/%dmax %dx%d, %lf*exp(%lf*i)",
							count_max, m, n, scale, rotation);
	//write sample image file extension based on parameters to generate it
	img_gen(m, n, count_max, scale, rotation, c1, c2, output_filename);
	//call to julia_gen.c to generate a Julia set sample image
	printf("\n");

	samples = grab_samples("./Samples");
	/*parse ./Samples/ for sample images meeting the naming convention...
	able to handle samples from previous runs for multiple comparisons*/

	for (i = 0; i < directories->num_dir; i++) {//loop trained images
		printf("Correlation with %s:\n", 
				trained_images[i].file_ext);//display trained image name
		for (j = 0; j < num_samples; j++) {//loop sample images
			if (trained_images[i].m == samples[j].m && trained_images[i].n 
				== samples[j].n) { //i.e. if dimensions match
				correlation = compare(trained_images[i], samples[j]);
				//computes correlation coefficient between trained and sample
				printf("\t%s: %lf\n", samples[j].file_ext, correlation);
				//display sample name and correlation to ith trained image
			}
		}
		printf("\n");
	}
	free(trained_images);
	free(samples);
	free(output_dir);
	for (i = 0; i < directories->num_dir; i++) {
		free(directories->images[i]);
	}
	free(directories->dir_ext);
  	free(directories->min);
  	free(directories->max);
  	free(directories->images);
  	free(directories->num_image);
  	free(directories->scale_flags);
  	free(directories);
  	free(trained_images);
  	//freeing malloc'ed arrays
	return 0;
}//end