#include "utils.h"
    

// DumpFieldToFile: dumps array into a file using RFS format


void DumpFieldToFile(int sx, int sy, int sz,
		     int ixStart, int ixEnd,
		     int iyStart, int iyEnd,
		     int izStart, int izEnd,
		     float d1, float d2, float d3,
		     float *arrP, char *fName) {

//PPL  int ix, iy, iz;
  int iy, iz;
  
  // header and binary file names
  
  char fNameHeader[128];
  strcpy(fNameHeader,fName);
  strcat(fNameHeader,".rfs");
  char fNameBinary[128];
  strcpy(fNameBinary,fNameHeader);
  strcat(fNameBinary,"@");

  // create header file in rfs format
  
  FILE *fp;
  fp=fopen(fNameHeader, "w+");
  fprintf(fp,"in=\"./%s\"\n", fNameBinary);
  fprintf(fp,"data_format=\"native_float\"\n");
  fprintf(fp,"esize=%lu\n", sizeof(float)); 
  fprintf(fp,"n1=%d\n",ixEnd-ixStart+1);
  fprintf(fp,"d1=%f\n",d1);
  fprintf(fp,"n2=%d\n",iyEnd-iyStart+1);
  fprintf(fp,"d2=%f\n",d2);
  fprintf(fp,"n3=%d\n",izEnd-izStart+1);
  fprintf(fp,"d3=%f\n",d3);
  fclose(fp);

  // create binary file
  
  fp=fopen(fNameBinary, "w+b");
  for (iz=izStart; iz<=izEnd; iz++)
    for (iy=iyStart; iy<=iyEnd; iy++) 
      fwrite((void *) (arrP +ind(ixStart,iy,iz)), sizeof(float), ixEnd-ixStart+1, fp);
  fclose(fp);
  int fsize=(ixEnd-ixStart+1)*(iyEnd-iyStart+1)*(izEnd-izStart+1);
  printf("created header (%s) and binary (%s) files\n",
	 fNameHeader,fNameBinary);
  printf("file has %d floats and %lu Mbytes from (%d,%d,%d) to (%d,%d,%d)\n",
	 fsize,fsize*sizeof(float)/(1024*1024),ixStart,iyStart,izStart,ixEnd,iyEnd,izEnd);
}


// DumpSlicePtr: prints the components of an SlicePtr variable


void DumpSlicePtr(SlicePtr p){
  printf("File %s contains time evolution of (%d:%d,%d:%d,%d:%d)\n",
	 p->fNameHeader,
	 p->ixStart, p->ixEnd, 
	 p->iyStart, p->iyEnd, 
	 p->izStart, p->izEnd);
}


// OpenSliceFile: open file in RFS format that will be continuously appended


SlicePtr OpenSliceFile(int ixStart, int ixEnd,
		       int iyStart, int iyEnd,
		       int izStart, int izEnd,
		       float dx, float dy, float dz, float dt,
		       char *fName) {
//PPL  char procName[128]="**(OpenSliceFile)**";
  SlicePtr ret;
  ret = (SlicePtr) malloc(sizeof(Slice));

  // verify slice direction

  if (ixStart==ixEnd)
    ret->direction=XSLICE;
  else if (iyStart==iyEnd)
    ret->direction=YSLICE;
  else if (izStart==izEnd)
    ret->direction=ZSLICE;
  else {
    ret->direction=FULL;
  }

  // header and binary file names
  
  strcpy(ret->fNameHeader,fName);
  strcat(ret->fNameHeader,".rfs");
  strcpy(ret->fNameBinary,FNAMEBINARYPATH);
  strcat(ret->fNameBinary,ret->fNameHeader);
  strcat(ret->fNameBinary,"@");

  // create header and binary files in rfs format
  
  ret->fpHead=fopen(ret->fNameHeader, "w+");
  ret->fpBinary=fopen(ret->fNameBinary, "w+");
  ret->ixStart=ixStart;
  ret->ixEnd=ixEnd;
  ret->iyStart=iyStart;
  ret->iyEnd=iyEnd;
  ret->izStart=izStart;
  ret->izEnd=izEnd;
  ret->itCnt=0;
  ret->dx=dx;
  ret->dy=dy;
  ret->dz=dz;
  ret->dt=dt;

  char sName[16];
  switch(ret->direction) {
  case XSLICE:
    strcpy(sName,"XSlice");
    break;
  case YSLICE:
    strcpy(sName,"YSlice");
    break;
  case ZSLICE:
    strcpy(sName,"ZSlice");
    break;
  case FULL:
    strcpy(sName,"Grid Section");
    break;
  }
  return(ret);
}

SlicePtr OpenSliceFile2(int ixStart, int ixEnd,
           int iyStart, int iyEnd,
           int izStart, int izEnd,
           float dx, float dy, float dz, float dt,
           char *fName) {
//PPL  char procName[128]="**(OpenSliceFile)**";
  SlicePtr ret;
  ret = (SlicePtr) malloc(sizeof(Slice));

  // verify slice direction

  if (ixStart==ixEnd)
    ret->direction=XSLICE;
  else if (iyStart==iyEnd)
    ret->direction=YSLICE;
  else if (izStart==izEnd)
    ret->direction=ZSLICE;
  else {
    ret->direction=FULL;
  }

  // header and binary file names
  
  strcpy(ret->fNameHeader,fName);
  strcat(ret->fNameHeader,".rfs");
  strcpy(ret->fNameBinary,FNAMEBINARYPATH);
  strcat(ret->fNameBinary,ret->fNameHeader);
  strcat(ret->fNameBinary,"@");

  // create header and binary files in rfs format
  
  ret->fpHead=fopen(ret->fNameHeader, "r+");
  ret->fpBinary=fopen(ret->fNameBinary, "r+");
  ret->ixStart=ixStart;
  ret->ixEnd=ixEnd;
  ret->iyStart=iyStart;
  ret->iyEnd=iyEnd;
  ret->izStart=izStart;
  ret->izEnd=izEnd;
  ret->itCnt=0;
  ret->dx=dx;
  ret->dy=dy;
  ret->dz=dz;
  ret->dt=dt;

  char sName[16];
  switch(ret->direction) {
  case XSLICE:
    strcpy(sName,"XSlice");
    break;
  case YSLICE:
    strcpy(sName,"YSlice");
    break;
  case ZSLICE:
    strcpy(sName,"ZSlice");
    break;
  case FULL:
    strcpy(sName,"Grid Section");
    break;
  }
  return(ret);
}


// DumpSliceFile: appends one array to an opened RFS file 


void DumpSliceFile2(int sx, int sy, int sz, int it,
		   float *arrP, SlicePtr p) {

      if(fseek(p->fpBinary, it * sizeof(float)  * (p->izEnd - p->izStart + 1) * (p->iyEnd - p->iyStart + 1) * (p->ixEnd-p->ixStart+1), SEEK_SET) == 1){
        printf("error reading slice!\n");
        exit(1);
      }

      fread((void *) (arrP),
	     sizeof(float),
	     (p->izEnd - p->izStart + 1) * (p->iyEnd - p->iyStart + 1) * (p->ixEnd-p->ixStart+1),
	     p->fpBinary);
}

void DumpSliceFile(int sx, int sy, int sz,
       float *arrP, SlicePtr p) {

//PPL  int ix, iy, iz;
  int iy, iz;
  
  // dump section to binary file
  for (iz=p->izStart; iz<=p->izEnd; iz++)
    for (iy=p->iyStart; iy<=p->iyEnd; iy++) 
      fwrite((void *) (arrP+ind(p->ixStart,iy,iz)),
       sizeof(float),
       p->ixEnd-p->ixStart+1,
       p->fpBinary);
//fsync(fileno(p->fpBinary));
  // increase it count
  
  p->itCnt++;
}


// CloseSliceFile: close file in RFS format that has been continuously appended


void CloseSliceFile(SlicePtr p){

  fprintf(p->fpHead,"in=\"%s\"\n", p->fNameBinary);
  fprintf(p->fpHead,"data_format=\"native_float\"\n");
  fprintf(p->fpHead,"esize=%lu\n", sizeof(float)); 
  switch(p->direction) {
  case XSLICE:
    fprintf(p->fpHead,"n1=%d\n",p->iyEnd-p->iyStart+1);
    fprintf(p->fpHead,"n2=%d\n",p->izEnd-p->izStart+1);
    fprintf(p->fpHead,"n3=%d\n",p->itCnt);
    fprintf(p->fpHead,"d1=%f\n",p->dy);
    fprintf(p->fpHead,"d2=%f\n",p->dz);
    fprintf(p->fpHead,"d3=%f\n",p->dt);
    break;
  case YSLICE:
    fprintf(p->fpHead,"n1=%d\n",p->ixEnd-p->ixStart+1);
    fprintf(p->fpHead,"n2=%d\n",p->izEnd-p->izStart+1);
    fprintf(p->fpHead,"n3=%d\n",p->itCnt);
    fprintf(p->fpHead,"d1=%f\n",p->dx);
    fprintf(p->fpHead,"d2=%f\n",p->dz);
    fprintf(p->fpHead,"d3=%f\n",p->dt);
    break;
  case ZSLICE:
    fprintf(p->fpHead,"n1=%d\n",p->ixEnd-p->ixStart+1);
    fprintf(p->fpHead,"n2=%d\n",p->iyEnd-p->iyStart+1);
    fprintf(p->fpHead,"n3=%d\n",p->itCnt);
    fprintf(p->fpHead,"d1=%f\n",p->dx);
    fprintf(p->fpHead,"d2=%f\n",p->dy);
    fprintf(p->fpHead,"d3=%f\n",p->dt);
    break;
  case FULL:
    fprintf(p->fpHead,"n1=%d\n",p->ixEnd-p->ixStart+1);
    fprintf(p->fpHead,"n2=%d\n",p->iyEnd-p->iyStart+1);
    fprintf(p->fpHead,"n3=%d\n",p->izEnd-p->izStart+1);
    fprintf(p->fpHead,"n4=%d\n",p->itCnt);
    fprintf(p->fpHead,"d1=%f\n",p->dx);
    fprintf(p->fpHead,"d2=%f\n",p->dy);
    fprintf(p->fpHead,"d3=%f\n",p->dz);
    fprintf(p->fpHead,"d4=%f\n",p->dt);
    break;
  }
  fclose(p->fpHead);
  fclose(p->fpBinary);
}


void CloseSliceFile2(SlicePtr p){


  fclose(p->fpHead);
  fclose(p->fpBinary);
}

// DumpSliceSummary: prints info of one array 


void   DumpSliceSummary(int sx, int sy, int sz,
			SlicePtr p,
			float dt, int it, float *arrP) {

  int ix, iy, iz;
  float maxP, minP, valP;
  maxP=minP=arrP[ind(p->ixStart,p->iyStart,p->izStart)];
  for (iz=p->izStart; iz<=p->izEnd; iz++)
    for (iy=p->iyStart; iy<=p->iyEnd; iy++) 
      for (ix=p->ixStart; ix<=p->ixEnd; ix++) {
	valP=arrP[ind(ix,iy,iz)];
	maxP=fmaxf(maxP,valP);
	minP=fminf(minP,valP);
      }
  printf("Slice of %s at iteration; %5.5d; time(s); %f; max; %9.2e; min; %9.2e; src; %9.2e;\n",
	 p->fNameHeader, it, dt*(float)it, maxP, minP, Source(dt,it));
}
