import React from 'react';
import ReactDOM from 'react-dom';
import { Link } from 'react-router-dom';

export default class Output extends React.Component {

	constructor(props) {
		super(props);
		this.state = {};
		this.querydata=this.querydata.bind(this);
	};

	querydata(){

		if(this.state.profile_result_Url == undefined){
			fetch('api/profiling/profile/'+window.batchid,{method:'GET'})
			.then(response => response.json())
			.then(result => this.setState(state=>({profile_result_Url:result.file})));
		}else if(this.state.png_file == undefined && this.state.profile_result_Url != undefined){

			fetch('api/profiling/dendrogram/'+window.batchid,{method:'GET'})
			.then(response => response.json())
			.then(result => this.setState(state=>({png_file:result.png_file,pdf_file:result.pdf_file,
				svg_file:result.svg_file,emf_file:result.emf_file,newick_file:result.newick_file })));

		}else if(this.state.png_file != undefined && this.state.profile_result_Url != undefined){

			this.setState(state=>({All_done:1}));
			clearInterval(this.interval);
		}

	}

	componentDidMount(){
		this.interval = setInterval(this.querydata,30000);
	}




    render() {

    	if(this.state.All_done == undefined){

    		return(<img src="https://www.cluecon.com/theme/img/pleasewait.gif" alt="Please waiting..." />);
    	
    	}else{
    		return (
    				<div id="url">
    				<a href={this.state.profile_result_Url}>profile file</a>
    				<div>Dendrogram:</div>
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
