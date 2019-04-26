import React from 'react';
import ReactDOM from 'react-dom';
import { Link } from 'react-router-dom';
import { Prompt } from 'react-router';
import Button from '@material-ui/core/Button';
import { withStyle } from '@material-ui/core/styles';
import ReplyIcon from '@material-ui/icons/Reply';
import CircularProgress from '@material-ui/core/CircularProgress';
import download from 'downloadjs';

export default class Dendrogram_view extends React.Component {

	constructor(props) {
		super(props);
		this.state = {};
		this.query_dengrogram = this.query_dengrogram.bind(this);
	};

	query_dengrogram(){
		
		if(this.state.png_file == undefined){
			fetch('api/dendrogram/dendrogram/' + window.clusteringID, { method: 'GET'})
			.then(response => response.json())
			.then(result => this.setState(state => ({
				png_file: result.png_file, 
				pdf_file: result.pdf_file,
				svg_file: result.svg_file, 
				emf_file: result.emf_file, 
				newick_file: result.newick_file })));
			
		}else {
			clearInterval(this.interval);
		}

	}

	componentDidMount(){
		this.interval = setInterval(this.query_dengrogram, 1000);
	}

	getIdFile(){
        download(window.clusteringID,'clustering_JobID.txt',"text/tab-separated-values");
    }

    render() {

    	if(this.state.png_file == undefined){

    		return(
	    		<div>
	    			<Prompt 
                            when={true} 
                            message="Are you sure to leave now?"/>
					<br />
					<br />
					<br />
					<div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                        <font size="6"> ID : {window.clusteringID}</font>
                        &nbsp;&nbsp;&nbsp;&nbsp;
                        <Button variant="contained" color="default" onClick={this.getIdFile}>
                            Get ID
                        </Button>
                    </div>
                    <br />
					<br />
					<div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
						<font size="6"> Please hold on ... </font>
					</div>
					<br />
					<br />
					<br />
					<div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
						<CircularProgress size={175} />
	                </div>
					<br />
					<br />
					<br />
				</div>
		);
    	
    	}else{
    		return (
				<div>
					<Prompt 
                            when={true} 
                            message="You are leaving the page. Please save ID to get result. Are you sure to leave now?"/>
					<br />
					<br />
					<div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                        <font size="6"> ID : {window.clusteringID}</font>
                        &nbsp;&nbsp;&nbsp;&nbsp;
                        <Button variant="contained" color="default" onClick={this.getIdFile}>
                            Get ID
                        </Button>
                    </div>
                    <br />
					<br />
					<div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
						<img src={this.state.svg_file} />
					</div>
					<br />
					<div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
						<font>Download</font>
					</div>
					<br />
					<div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
						<a download href={this.state.png_file} style={{ textDecoration:'none' }}>
							<Button variant="contained" color="default">Png</Button>
						</a>
						&nbsp;&nbsp;&nbsp;&nbsp;
						<a download href={this.state.pdf_file} style={{ textDecoration:'none' }}>
							<Button variant="contained" color="default">Pdf</Button>
						</a>
						&nbsp;&nbsp;&nbsp;&nbsp;
						<a download href={this.state.svg_file} style={{ textDecoration:'none' }}>
							<Button variant="contained" color="default">Svg</Button>
						</a>
						&nbsp;&nbsp;&nbsp;&nbsp;
						<a download href={this.state.emf_file} style={{ textDecoration:'none' }}>
							<Button variant="contained" color="default">emf</Button>
						</a>
						&nbsp;&nbsp;&nbsp;&nbsp;
						<a download href={this.state.newick_file} style={{ textDecoration:'none' }}>
							<Button variant="contained" color="default">newick</Button>
						</a>
					</div>
					<br />
					<br />
					<br />
					<div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
						<Link to="/upload_profile" style={{ textDecoration:'none' }}>
							<Button variant="contained" color="default">
								<ReplyIcon />
								&nbsp;&nbsp;
								Back
							</Button>
						</Link>
					</div>
					<br />
					<br />
				</div>
			);
		}
    }
}
