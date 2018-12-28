import React from 'react';
import ReactDOM from 'react-dom';
import Navigation from './Navigation.jsx';
import { Link } from 'react-router-dom';
import Paper from '@material-ui/core/Paper';
import Tabs from '@material-ui/core/Tabs';
import Tab from '@material-ui/core/Tab';
import Button from '@material-ui/core/Button';
import ReplyIcon from '@material-ui/icons/Reply';
import DownloadIcon from '@material-ui/icons/CloudDownload';

export default class Example extends React.Component {

	constructor(props) {
		super(props);

        this.state = {};
		this.query_demo = this.query_demo.bind(this);
    }

	query_demo(){

        fetch('api/profiling/profile/816e94d3-2ad1-4886-95b1-5947c4a333a6/', { method:'GET'})
			.then(response => response.json())
			.then(result => this.setState(state => ({
                profile_result_all: result.file,
                profile_result_zip: result.zip })));

		fetch('api/dendrogram/dendrogram/816e94d3-2ad1-4886-95b1-5947c4a333a6/', { method: 'GET'})
            .then(response => response.json())
            .then(result => this.setState(state => ({
                png_file: result.png_file, 
                pdf_file: result.pdf_file,
                svg_file: result.svg_file, 
                emf_file: result.emf_file, 
                newick_file: result.newick_file })));
	}

	componentDidMount(){
		this.query_demo();
	}

    render() {
    	return (
			<div id="url">
				<Navigation value={3}/>
				<br />
				<br />
	            <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
	            	<a download href={this.state.profile_result_zip} 
                    style={{ textDecoration:'none' }}>
	                    <Button variant="contained" color="default">
		                    Download profiles (.zip)
		                    &nbsp;&nbsp;
		                    <DownloadIcon />
	                    </Button>
	                </a>
                    &nbsp;&nbsp;&nbsp;&nbsp;
                    <a download href={this.state.profile_result_all} 
                    style={{ textDecoration:'none' }}>
		                <Button variant="contained" color="default">
		                    Download profiles (.tsv)
		                    &nbsp;&nbsp;
		                    <DownloadIcon />
	                    </Button>
	                </a>
	            </div>
				<br />
				<br />
	            <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
	                <img src={require('./static/dendrogram_example.svg')} />
	            </div>
	            <br />
	            <br />
				<br />
	            <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
	            	<font size="4">Download </font> 
	            	&nbsp;&nbsp;&nbsp;&nbsp;
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
			</div>
		);

        
    }
}