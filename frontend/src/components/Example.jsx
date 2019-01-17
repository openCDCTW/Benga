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
    }

    render() {
    	return (
			<div id="url">
				<Navigation value={3}/>
				<br />
				<br />
	            <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
	            	<a download href={require('./static/example_result/profiles-example.zip')} 
                    style={{ textDecoration:'none' }}>
	                    <Button variant="contained" color="default">
		                    Download profiles (.zip)
		                    &nbsp;&nbsp;
		                    <DownloadIcon />
	                    </Button>
	                </a>
                    &nbsp;&nbsp;&nbsp;&nbsp;
                    <a download href={require('./static/example_result/profiles-example.tsv')} 
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
	                <img src={require('./static/example_result/dendrogram_example.svg')} />
	            </div>
	            <br />
	            <br />
				<br />
	            <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
	            	<font size="4">Download </font> 
	            	&nbsp;&nbsp;&nbsp;&nbsp;
	            	<a download href={require('./static/example_result/dendrogram_example.png')} 
	            	style={{ textDecoration:'none' }}>
	            		<Button variant="contained" color="default">Png</Button>
	            	</a>
	                &nbsp;&nbsp;&nbsp;&nbsp;
	                <a download href={require('./static/example_result/dendrogram_example.pdf')} 
	                style={{ textDecoration:'none' }}>
	                	<Button variant="contained" color="default">Pdf</Button>
	                </a>
	                &nbsp;&nbsp;&nbsp;&nbsp;
	                <a download href={require('./static/example_result/dendrogram_example.svg')} 
	                style={{ textDecoration:'none' }}>
	                	<Button variant="contained" color="default">Svg</Button>
	                </a>
	                &nbsp;&nbsp;&nbsp;&nbsp;
	                <a download href={require('./static/example_result/dendrogram_example.emf')} 
	                style={{ textDecoration:'none' }}>
	                	<Button variant="contained" color="default">emf</Button>
	                </a>
	                &nbsp;&nbsp;&nbsp;&nbsp;
	                <a download href={require('./static/example_result/dendrogram_example.newick')} 
	                style={{ textDecoration:'none' }}>
	                	<Button variant="contained" color="default">newick</Button>
	                </a>
	            </div>
	            <br />
	            <br />
			</div>
		);

        
    }
}