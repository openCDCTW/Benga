import React from 'react';
import ReactDOM from 'react-dom';

export default class Dendrogram extends React.Component {

	constructor(props) {
		super(props);
		this.state = {};
		this.query_dengrogram = this.query_dengrogram.bind(this);
	};

	query_dengrogram(){

		if(this.state.png_file == undefined){
			fetch('api/profiling/dendrogram/'+window.batchid,{method:'GET'})
			.then(response => response.json())
			.then(result => this.setState(state=>({png_file:result.png_file, pdf_file:result.pdf_file,
				svg_file:result.svg_file, emf_file:result.emf_file, newick_file:result.newick_file })));
			
		}else {
			clearInterval(this.interval);
		}

	}

	componentDidMount(){
		this.query_dengrogram();
		this.interval = setInterval(this.query_dengrogram,30000);
	}



    render() {

    	if(this.state.png_file == undefined){

    		return(<img src="http://i.imgur.com/1J9vdDD.gif" alt="Please waiting..." />);
    	
    	}else{
    		return (
    				<div id="url">
    				<img src={this.state.svg_file} />
    				<div>Download Dendrogram:</div>
    				<a download href={this.state.png_file}>png</a><br />
    				<a download href={this.state.pdf_file}>pdf</a><br />
    				<a download href={this.state.svg_file}>svg</a><br />
    				<a download href={this.state.emf_file}>emf</a><br />
    				<a download href={this.state.newick_file}>newick</a><br />
    				</div>
        	);
    	}
        
    }
}
