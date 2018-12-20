import React from 'react';
import ReactDOM from 'react-dom';
import Navigation from './Navigation.jsx';
import { Link } from 'react-router-dom';
import Paper from '@material-ui/core/Paper';
import Tabs from '@material-ui/core/Tabs';
import Tab from '@material-ui/core/Tab';
import Button from '@material-ui/core/Button';
import { withStyle } from '@material-ui/core/styles';
import ReplyIcon from '@material-ui/icons/Reply';
import DownloadIcon from '@material-ui/icons/CloudDownload';


export default class Example extends React.Component {

	constructor(props) {
		super(props);

        this.state = {};
		this.query_demo = this.query_demo.bind(this);
    }

	query_demo(){
		fetch('api/dendrogram/dendrogram/c327bb9e-2dea-4fe3-baee-88d0bcf4d93a/', { method: 'GET'})
            .then(response => response.json())
            .then(result => this.setState(state => ({
                svg_file: result.svg_file })));
	}

	componentDidMount(){
		this.query_demo();
	}

    render() {

    	return (
			<div id="url">
				<Navigation value={2}/>
				<br />
				<br />
	            <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                    <Button variant="contained" color="default">
	                    Download profiles (.zip)
	                    &nbsp;&nbsp;
	                    <DownloadIcon />
                    </Button>
                    &nbsp;&nbsp;&nbsp;&nbsp;
	                <Button variant="contained" color="default">
	                    Download profiles (.tsv)
	                    &nbsp;&nbsp;
	                    <DownloadIcon />
                    </Button>
	            </div>
				<br />
				<br />
	            <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
	                <img src={this.state.svg_file} />
	            </div>
	            <br />
	            <br />
				<br />
	            <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
	            	<font size="4">Download </font> 
	            	&nbsp;&nbsp;&nbsp;&nbsp;
	            	<Button variant="contained" color="default">Png</Button>
	                &nbsp;&nbsp;&nbsp;&nbsp;
	                <Button variant="contained" color="default">Pdf</Button>
	                &nbsp;&nbsp;&nbsp;&nbsp;
	                <Button variant="contained" color="default">Svg</Button>
	                &nbsp;&nbsp;&nbsp;&nbsp;
	                <Button variant="contained" color="default">emf</Button>
	                &nbsp;&nbsp;&nbsp;&nbsp;
	                <Button variant="contained" color="default">newick</Button>
	            </div>
	            <br />
	            <br />
			</div>
		);

        
    }
}
