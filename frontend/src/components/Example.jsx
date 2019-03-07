import React from 'react';
import ReactDOM from 'react-dom';
import { Link } from 'react-router-dom';
import Paper from '@material-ui/core/Paper';
import Tabs from '@material-ui/core/Tabs';
import Tab from '@material-ui/core/Tab';
import Button from '@material-ui/core/Button';
import ReplyIcon from '@material-ui/icons/Reply';
import DownloadIcon from '@material-ui/icons/CloudDownload';
import download from 'downloadjs';

export default class Example extends React.Component {

    download_zip(token){
		fetch(require('./static/Example_result/profiles-example.zip'),{ 
			method:'GET',
			headers: {
				'Authorization':token
			}})
		.then(response => response.blob())
		.then(blob => {
			download(blob,'profiles-example');
		});
	}

    download_tsv(token){
		fetch(require('./static/Example_result/profiles-example.tsv'),{ 
			method:'GET',
			headers: {
				'Authorization':token
			}})
		.then(response => response.blob())
		.then(blob => {
			download(blob,'profiles-example.tsv',"text/tab-separated-values");
		});
	}

	download_png(token){
		fetch(require('./static/Example_result/dendrogram_example.png'),{ 
			method:'GET',
			headers: {
				'Authorization':token
			}})
		.then(response => response.blob())
		.then(blob => {
			download(blob,'dendrogram_example');
		});
	}

	download_pdf(token){
		fetch(require('./static/Example_result/dendrogram_example.pdf'),{ 
			method:'GET',
			headers: {
				'Authorization':token
			}})
		.then(response => response.blob())
		.then(blob => {
			download(blob,'dendrogram_example');
		});
	}

	download_svg(token){
		fetch(require('./static/Example_result/dendrogram_example.svg'),{ 
			method:'GET',
			headers: {
				'Authorization':token
			}})
		.then(response => response.blob())
		.then(blob => {
			download(blob,'dendrogram_example');
		});
	}

	download_emf(token){
		fetch(require('./static/Example_result/dendrogram_example.emf'),{ 
			method:'GET',
			headers: {
				'Authorization':token
			}})
		.then(response => response.blob())
		.then(blob => {
			download(blob,'dendrogram_example.emf',"image/emf");
		});
	}

	download_newick(token){
		fetch(require('./static/Example_result/dendrogram_example.newick'),{ 
			method:'GET',
			headers: {
				'Authorization':token
			}})
		.then(response => response.blob())
		.then(blob => {
			download(blob,'dendrogram_example.newick');
		});
	}

    render() {
    	return (
			<div id="url">
				<br />
				<br />
	            <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                    <Button variant="contained" color="default" onClick={this.download_zip.bind(this)}>
	                    Download profiles (.zip)
	                    &nbsp;&nbsp;
	                    <DownloadIcon />
                    </Button>
                    &nbsp;&nbsp;&nbsp;&nbsp;
	                <Button variant="contained" color="default" onClick={this.download_tsv.bind(this)}>
	                    Download profiles (.tsv)
	                    &nbsp;&nbsp;
	                    <DownloadIcon />
                    </Button>
	            </div>
				<br />
				<br />
	            <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
	                <img src={require('./static/Example_result/dendrogram_example.svg')} />
	            </div>
	            <br />
	            <br />
				<br />
	            <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
	            	<font size="4">Download </font> 
	            	&nbsp;&nbsp;&nbsp;&nbsp;
	            	<Button variant="contained" color="default" onClick={this.download_pdf.bind(this)}>
	                	Pdf
	                </Button>
	                &nbsp;&nbsp;&nbsp;&nbsp;
	                <Button variant="contained" color="default" onClick={this.download_emf.bind(this)}>
	                	emf
	                </Button>
	                &nbsp;&nbsp;&nbsp;&nbsp;
	                <Button variant="contained" color="default" onClick={this.download_svg.bind(this)}>
	               		Svg
	               	</Button>
	                &nbsp;&nbsp;&nbsp;&nbsp;
	            	<Button variant="contained" color="default" onClick={this.download_png.bind(this)}>
	            		Png
	            	</Button>
	                &nbsp;&nbsp;&nbsp;&nbsp;
	                <Button variant="contained" color="default" onClick={this.download_newick.bind(this)}>
	                	newick
	                </Button>
	            </div>
	            <br />
	            <br />
			</div>
		);

        
    }
}