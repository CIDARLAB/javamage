package utils;
import java.util.ArrayList;
import java.util.List;

import com.panayotis.gnuplot.JavaPlot;
import com.panayotis.gnuplot.layout.AutoGraphLayout;
import com.panayotis.gnuplot.plot.AbstractPlot;
import com.panayotis.gnuplot.plot.DataSetPlot;
import com.panayotis.gnuplot.style.Style;


public class Plot {
	
	private JavaPlot	jplot;
	private ArrayList<SubGraph> graphs;

	public Plot(Double[] Ydata, String title) {
		this();
		this.addGraph(Ydata, title);
	}
	
	public	Plot(){
		// Initialize
		jplot = new JavaPlot();
		graphs = new ArrayList<SubGraph>();
		
		// Set to a column style
		AutoGraphLayout layout =((AutoGraphLayout) jplot.getLayout());
		layout.setColumns(1);
		//layout.setRows(3);
		
		// Set aesthetic properties
		jplot.set("key", "outside right spacing 0.75");
		jplot.set("screen-coordinate", "min 0,0 max 1,1");
		jplot.set("border", "3");
		jplot.set("xtics", "nomirror");
		jplot.set("ytics", "nomirror");
		jplot.set("bmargin", "1");
		jplot.set("tmargin", "2");
	}
	
	public void  addGraph(List<Double[]> Ydatas, List<String> titles){

		ArrayList<DataSetPlot> data = new ArrayList<DataSetPlot>();
		
		for (int ii = 0 ; ii< Ydatas.size() ;ii++){

			Double[] yy = Ydatas.get(ii);

			double[][] dataArray = new double[yy.length][2];

			// Populate the dataArray
			for( int jj = 0; jj <yy.length; jj++ ) {
				dataArray[jj][0]= (double)jj; 
				dataArray[jj][1]= yy[jj];
			}

			DataSetPlot set = new DataSetPlot(dataArray);
			set.setTitle(titles.get(ii));
			data.add(set);
		}
		
		SubGraph sg = new SubGraph(data); 
		graphs.add(sg);

	}

	public void addGraph(List<Double[]> Xdatas, List<Double[]> Ydatas, List<String> titles) {

		ArrayList<DataSetPlot> data = new ArrayList<DataSetPlot>();

		for (int ii = 0 ; ii< Ydatas.size() ;ii++){

			Double[] xx = Xdatas.get(ii);
			Double[] yy = Ydatas.get(ii);

			double[][] dataArray = new double[xx.length][2];

			// Populate the dataArray
			for( int jj = 0; jj <yy.length; jj++ ) {
				dataArray[jj][0]= xx[jj]; 
				dataArray[jj][1]= yy[jj];
			}

			DataSetPlot set = new DataSetPlot(dataArray);
			set.setTitle(titles.get(ii));
			data.add(set);
		}
		
		SubGraph sg = new SubGraph(data); 
		graphs.add(sg);

//		for (DataSetPlot set: data) {
//			Graph gr =  new Graph();
//			gr.add(set);
//			graphs.add(gr);
//		}

	}

	public void addGraph(Double[] Xdata, Double[] Ydata, String title) {


		double[][] dataArray = new double[Ydata.length][2];

		for( int ii = 0; ii <Ydata.length; ii++ ) {
			dataArray[ii][0]= Xdata[ii]; 
			dataArray[ii][1]= Ydata[ii];
		}

		DataSetPlot set = new DataSetPlot(dataArray);
		set.setTitle(title);
		ArrayList<DataSetPlot> data = new ArrayList<DataSetPlot>();
		data.add(set);
		
		SubGraph sg = new SubGraph(data); 
		graphs.add(sg);
	}

	public void addGraph(Double[] Ydata, String title) {


		// Populate the dataArray
		double[][] dataArray	= new double[Ydata.length][2];
		for( int ii = 0; ii <Ydata.length; ii++ ) {
			dataArray[ii][0]= (double)ii; 
			dataArray[ii][1]= Ydata[ii];
		}

		DataSetPlot set = new DataSetPlot(dataArray);
		set.setTitle(title);
		
		ArrayList<DataSetPlot> data = new ArrayList<DataSetPlot>();
		data.add(set);
		SubGraph sg = new SubGraph(data); 
		graphs.add(sg);
	}

	public void title(String titleString) {
		jplot.setMultiTitle(titleString);
	}

	public void setToLines(){
		for (SubGraph sg : graphs ) {
			for (DataSetPlot set : sg.getDataSets()) {
				((AbstractPlot) set).getPlotStyle().setStyle(Style.LINES);
			}
		}
	}

	public void setToPoints(){
		for (SubGraph sg : graphs ) {
			for (DataSetPlot set : sg.getDataSets()) {
				((AbstractPlot) set).getPlotStyle().setStyle(Style.POINTS);
			}
		}
	}

	public void setToLinePoints(){
		for (SubGraph sg : graphs ) {
			for (DataSetPlot set : sg.getDataSets()) {
				((AbstractPlot) set).getPlotStyle().setStyle(Style.LINESPOINTS);
			}
		}
	}

	public void xtitle(String titleString) {
		jplot.getAxis("x").setLabel(titleString);
	}

	public void ytitle(String titleString) {
		jplot.getAxis("y").setLabel(titleString);
	}

	public void draw() {		
		
		boolean first=true;
		for ( SubGraph sg: graphs ){
			if (!first) { jplot.newGraph(); } else {first= false;}	
			this.addSubGraph(sg);	
		}
		jplot.plot();
	}

	private void addSubGraph(SubGraph sg){
		for (DataSetPlot dsp : sg.getDataSets()) {
			jplot.addPlot(dsp);
		}
	}
	
	public static void quickPlot(Double[] dd,String title){
		Plot plotwindow = new Plot(dd,title);
		plotwindow.setToLinePoints();
		plotwindow.draw();

	}
	
	public class SubGraph{
		
		private List<DataSetPlot>	data;
		
		public SubGraph(List<DataSetPlot> ll){
			
			this.data = ll;
		}
		
		public List<DataSetPlot> getDataSets(){
			
			return data;
		}
	}


}
