**	Plot distribution of MRS by group;

data temp;
	Dose = 'Placebo';
		MRS = 0; Frequency = .09; output;
		MRS = 1; Frequency = .13; output;
		MRS = 2; Frequency = .09; output;
		MRS = 3; Frequency = .13; output;
		MRS = 4; Frequency = .16; output;
		MRS = 5; Frequency = .09; output;
		MRS = 6; Frequency = .31; output;
	Dose = 'Dose 1';
		MRS = 0; Frequency = .12; output;
		MRS = 1; Frequency = .16; output;
		MRS = 2; Frequency = .11; output;
		MRS = 3; Frequency = .15; output;
		MRS = 4; Frequency = .15; output;
		MRS = 5; Frequency = .12; output;
		MRS = 6; Frequency = .19; output;
	Dose = 'Dose 2';
		MRS = 0; Frequency = .15; output;
		MRS = 1; Frequency = .20; output;
		MRS = 2; Frequency = .10; output;
		MRS = 3; Frequency = .16; output;
		MRS = 4; Frequency = .15; output;
		MRS = 5; Frequency = .10; output;
		MRS = 6; Frequency = .14; output;
	Dose = 'Dose 3';
		MRS = 0; Frequency = .17; output;
		MRS = 1; Frequency = .24; output;
		MRS = 2; Frequency = .15; output;
		MRS = 3; Frequency = .09; output;
		MRS = 4; Frequency = .20; output;
		MRS = 5; Frequency = .04; output;
		MRS = 6; Frequency = .11; output;
run;

ods pdf file = "C:\Users\rreeder\Desktop\AdaptiveDoseFinding\MRS distribution by dose.pdf";
	proc template;
	define statgraph sgdesign;
	dynamic _MRS _FREQUENCY _DOSE;
	begingraph;
	   entrytitle halign=center 'Distribution of MRS by Dosing Group';
	   layout lattice / rowdatarange=data columndatarange=data rowgutter=10 columngutter=10;
	      layout overlay / xaxisopts=( discreteopts=( tickvaluefitpolicy=splitrotate)) yaxisopts=( label=('Probability'));
	         barchart category=_MRS response=_FREQUENCY / group=_DOSE name='bar' stat=mean groupdisplay=Cluster clusterwidth=0.85;
	         discretelegend 'bar' / opaque=false border=true halign=left valign=bottom displayclipped=true across=1 order=rowmajor location=inside autoalign=(topright topleft bottomright bottomleft top bottom right left);
	      endlayout;
	   endlayout;
	endgraph;
	end;
	run;

	proc sgrender data=WORK.TEMP template=sgdesign;
	dynamic _MRS="MRS" _FREQUENCY="FREQUENCY" _DOSE="DOSE";
	run;
ods pdf close;