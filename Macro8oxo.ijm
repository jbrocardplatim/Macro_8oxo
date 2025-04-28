//Macro used to measure intensity of 8oxo labelling in tom+ cells from brain slices
//Needs the Cellpose wrapper installed : https://github.com/BIOP/ijl-utilities-wrappers
//@Jacques Brocard for Juan Ren & Sabine Richard (IGFL), 2025

//---INITIALIZATION
run("Set Measurements...", "area mean redirect=None decimal=3");
cell_radius=15; //Mean cell radius in pixels
nb_z=5; //Nb of planes
n_c1="8oxo"; //Labelling to be detected in 1st channel
n_c2="tom"; //Labelling to be detected in 2nd channel
cellpose_path="D:\\DATA-USERS\\PLATIM\\JAC\\cellposenv";

//---PICK DIRECTORY and analyzes each .czi stack
dir = getDirectory("Choose the Directory");
list = getFileList(dir);

for (i=0; i<list.length; i++){
    ext=substring(list[i],lengthOf(list[i])-4,lengthOf(list[i]));
    if (ext==".czi"){
    	print(list[i]);
    	run("Bio-Formats Windowless Importer", "open=["+dir+list[i]+"]");
    	getDimensions(w, h, channels, slices, frames);
    	nb_c=channels;
    	roiManager("Reset");
    	main(list[i]);
    	close();
    }
    selectWindow("Log");
	saveAs("Text",dir+"Results_8oxo_v8.txt");
} 


function main(t){

	//--- INITIALIZE STACK
	mean=newArray(nb_z);
	t=getTitle();
	dir=getDirectory("image");
		
	stack=substring(t,0,lengthOf(t)-4)+"_stack.tif";
	temp=substring(t,0,lengthOf(t)-4)+"_temp.tif";
	
	//Automatic projection of the three most intense 8oxo z planes
	for (s=0;s<floor(nSlices/nb_c);s++){
		setSlice(1+s*nb_c);
		getStatistics(area, mean[s]);
	}
	maxi=0;
	for (s=1;s<nb_z;s++){
		if (mean[s]>mean[maxi]) maxi=s;
	}
	deb=maxi-1;
	fin=maxi+1;
	if (deb<0) {
		deb++; fin++;
	}
	if (fin>nb_z) {
		deb--; fin--;
	}
	run("Z Project...", "start="+deb+" stop="+fin+" projection=[Max Intensity]");
	
	//Save reference original stack	
	rename(stack);
	saveAs("Tiff",dir+stack);
	
	//Smooth and subtract background before image analysis
	//Also, saturate 0.1% of pixels for intensity normalization
	setSlice(1);
	run("Select All");
	run("Properties...", "channels=2 slices=1 frames=1 pixel_width=1 pixel_height=1 voxel_depth=1 unit=pixel");
	run("Smooth","stack");
	run("Subtract Background...", "rolling="+2*cell_radius+" stack");
	saturate_slice(1,0.001);
	saturate_slice(2,0.001);
	
	//Save analysis-ready stack as a "_temp.tif" file
	roi_file=dir+substring(stack,0,lengthOf(stack)-4)+"_8oxo.zip";
	saveAs("Tiff",dir+temp);
	
	//Automatic detection of cell contours using Cellpose with channels 1 & 2 labellings
	//and extract ROIs of automatically detected cells (extract_cellpose) and more (add_roi)
	run("Cellpose ...", "env_path="+cellpose_path+" env_type=venv model=cyto3 model_path=path\\to\\own_cellpose_model diameter="+2*cell_radius+" ch1=1 ch2=2 additional_flags=");
	rename("tom");
	extract_cellpose("tom");
	add_rois(temp);
	roiManager("Deselect");
	roiManager("Save",roi_file);
	
	selectWindow(temp);
	//Classify cells based on individual labelling intensities
	classify_cellpose();
	close(); 
}

function add_rois(te){
	//Adds ROIs for brightly labelled tom+ cells that may have gone undetected
	selectWindow(te);
	setSlice(2);
	run("Duplicate...", "title=2 channels=2");
	setForegroundColor(0,0,0);
	selectWindow("ROI Manager");
	roiManager("Fill");
	run("Enhance Contrast", "saturated=0.5");
	getMinAndMax(min, max);
	setThreshold(max, 65535);
	setOption("BlackBackground", true);
	run("Convert to Mask");
	run("Options...", "iterations=2 count=1 black do=Dilate");
	run("Fill Holes");
	run("Watershed");
	cell_size=PI*cell_radius*cell_radius;
	run("Analyze Particles...", "size="+cell_size+"-Infinity circularity=0.50-1.00 exclude add");	
	close("2");
}

function extract_cellpose(te){
	//Extract ROIs of proper size from mask of automatically detected cells
	selectWindow(te);
	setAutoThreshold("Default dark");
	run("Select All");
	getStatistics(area, mean, min, max);
	cell_size=PI*cell_radius*cell_radius/2;
	for (i=1;i<max;i=i+0.99){
		setThreshold(i,i+1);
		run("Analyze Particles...", "size="+cell_size+"-Infinity exclude add");	
	}
	close(te);
}


function saturate_slice(s,percent) { 
	//Modify image display to indicated percentage of saturated pixels, before applying LUT
	setSlice(s);
	run("Select All");
	getHistogram(values, counts, 256);
	
	beg=0; area=0;
	while (beg<=255) {
		area=area+counts[beg];
		beg++;
	}
	
	beg=1;
	while(counts[beg+1]>counts[beg]) beg++;
	end=255;
	while(sum_counts(end)/area<percent) end--;
	temp_end=values[end];
	if (temp_end<15000) temp_end=15000;
	
	setMinAndMax(values[beg],temp_end);
	run("Apply LUT", "slice");
	
}


function sum_counts(e){
	sum=0;
	for (f=255;f>=e;f--){
		sum=sum+counts[f];
	}
	return sum;
}


function classify_cellpose(){

	//At this step, the ROI Manager contains all detected cell contours
	c1=0; c2=0;
	th_c1=12000; 
	th_c2=14000;  
	print("ROI\t Area\t Mean");
	
	nROIs=roiManager("Count"); pos=1;
	area=newArray(3); mean=newArray(3); 
	//For each ROI = cell contour
	for (r=0;r<nROIs;r++){
		roiManager("Select",0);
		for (s=1;s<=2;s++){
			setSlice(s);
			//Area and Mean in each channel is measured...
			getStatistics(area[s], mean[s]);
		}
		col="black"; c1pos=false; c2pos=false; 

		//... and means are compared to thresholds to establish identity
		if (mean[1]>th_c1) c1pos=true; //8oxo-positive cell
		if (mean[2]>th_c2) c2pos=true; //tom-positive cell

		if (!c1pos & c2pos) {
			c2++; col="red"; //tom+ 8oxo- cells will be red
		}
		if (c1pos & c2pos) {
			c1++; col="green"; //tom+ 8oxo+ cells will be green
			print(pos+"\t "+area[1]+"\t "+floor(mean[1]));
		}

		roiManager("Rename", pos+"-"+col);
		Roi.setStrokeColor(col);
		Roi.setPosition("none");
		if (col!="black") {
			roiManager("Add");
			pos++;
		}
		roiManager("Delete");
	}
	
	//Save color-coded ROI Manager and print % tom+ 8oxo+ cells
	pos=pos-1;
	if (pos>=0){
		nROIs=roiManager("Count");
		for (r=0;r<nROIs;r++){
			roiManager("Select",r);
			Roi.setPosition("none");
		}
		roiManager("Save",roi_file);
		print(t+"\t "+pos+"\t "+floor(1000*c1/pos)/10+"\t % 8oxo+ tom-pos cells");
	}
}