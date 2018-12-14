import React from 'react';
import ReactDOM from 'react-dom';
import { Link } from 'react-router-dom';
import Paper from '@material-ui/core/Paper';
import Tabs from '@material-ui/core/Tabs';
import Tab from '@material-ui/core/Tab';
import Button from '@material-ui/core/Button';
import { withStyle } from '@material-ui/core/styles';
import ReplyIcon from '@material-ui/icons/Reply';
import DownloadIcon from '@material-ui/icons/CloudDownload';


export default class QueryData extends React.Component {

	constructor(props) {
		super(props);

        this.state = {};
		this.query_profile = this.query_profile.bind(this);
        this.query_dendrogram = this.query_dendrogram.bind(this);
    }

	query_profile(){

        fetch('api/profiling/profile/' + window.queryID, { method:'GET'})
            .then(response => response.json())
            .then(result => this.setState(state => ({
                profile_result_all: result.file,
                profile_result_zip: result.zip,})));
	}

    query_dendrogram(){

        fetch('api/dendrogram/dendrogram/' + window.queryID, { method: 'GET'})
            .then(response => response.json())
            .then(result => this.setState(state => ({
                png_file: result.png_file, 
                pdf_file: result.pdf_file,
                svg_file: result.svg_file, 
                emf_file: result.emf_file, 
                newick_file: result.newick_file, })));
    }

    cleanID(){
        window.queryID = "";
    }

	componentDidMount(){
		this.query_profile();
        this.query_dendrogram();
	}

    render() {

    	if(this.state.profile_result_all == undefined){

    		return(
    			<div>
                    <Paper square>
                        <Tabs value={false} centered>
                            <Tab label=" " disabled/>
                        </Tabs>
                    </Paper>
                    <br />
                    <br />
                    <br />
                    <br />
                    <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                        <font size="4"> Data not found. </font>
                    </div>
                    <br />
                    <br />
                    <br />
                    <br />
                    <div style={{ display:'flex', justifyContent:'center', alignItems:'center' }}>
                        <Link to="/" style={{ textDecoration:'none' }}>
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
    	
    	}else{
    		return (
    				<div id="url">
                        <Paper square>
                            <Tabs centered>
                                <Tab label=" "/>
                            </Tabs>
                        </Paper>
    					<br />
                        <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                            <a download href={this.state.profile_result_all} 
                             style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                Download profiles (.tsv)
                                &nbsp;&nbsp;
                                <DownloadIcon />
                                </Button>
                            </a>
                            &nbsp;&nbsp;&nbsp;&nbsp;
                            <a download href={this.state.profile_result_zip} 
                             style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                Download profiles (.zip)
                                &nbsp;&nbsp;
                                <DownloadIcon />
                                </Button>
                            </a>
                        </div>
    					<br />
    					<br />
                        <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                            <img src={this.state.svg_file} />
                        </div>
                        <br />
                        <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                            <font size="4">Download Dendrogram</font>
                        </div>
                        <br />
                        <br />
                        <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                            
                            <a download href={this.state.png_file} style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                Png 
                                </Button>
                            </a>
                            &nbsp;&nbsp;&nbsp;&nbsp;
                            <a download href={this.state.pdf_file} style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                Pdf
                                </Button>
                            </a>
                            &nbsp;&nbsp;&nbsp;&nbsp;
                            <a download href={this.state.svg_file} style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                Svg
                                </Button>
                            </a>
                            &nbsp;&nbsp;&nbsp;&nbsp;
                            <a download href={this.state.emf_file} style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                emf
                                </Button>
                            </a>
                            &nbsp;&nbsp;&nbsp;&nbsp;
                            <a download href={this.state.newick_file} style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                newick
                                </Button>
                            </a>
                        </div>
                        <br />
                        <br />
                        <div style={{ display:'flex', justifyContent:'center', alignItems:'center' }}>
                            <Link to="/" style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default" onClick={this.cleanID.bind(this)}>
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
