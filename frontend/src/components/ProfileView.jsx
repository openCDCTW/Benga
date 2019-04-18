import React from 'react';
import ReactDOM from 'react-dom';
import { Link } from 'react-router-dom';
import Button from '@material-ui/core/Button';
import { withStyle } from '@material-ui/core/styles';
import ReplyIcon from '@material-ui/icons/Reply';
import DownloadIcon from '@material-ui/icons/CloudDownload';
import CircularProgress from '@material-ui/core/CircularProgress';
import download from 'downloadjs';
import { Prompt } from 'react-router';

export default class Profile_view extends React.Component {

	constructor(props) {
		super(props);

        let fileName_str="", i=0, j=0;
        for (i; i < window.fileName.length; i++){
            if(i < window.fileName.length-1){
                fileName_str += window.fileName[i] + ", "
            }else{
                fileName_str += window.fileName[i]
            }
        };

        this.state = { fileName :fileName_str };
		this.query_result = this.query_result.bind(this);
    }

	query_result(){

		if(this.state.profile_result_zip == undefined){
			fetch('api/profiling/profile/' + window.batchid, { method:'GET'})
			.then(response => response.json())
			.then(result => this.setState(state => ({
                profile_result_zip: result.zip })));
		}else{
			clearInterval(this.interval);
		}

	}

    getIdFile(){
        download(window.batchid,'BatchId.txt',"text/tab-separated-values");
    }

	componentDidMount(){
		this.interval = setInterval(this.query_result, 60000);
	}

    render() {

    	if(this.state.profile_result_zip == undefined){

    		return(
    			<div>
                    <Prompt 
                        when={true} 
                        message="You are leaving the page. Please save ID to get result. Are you sure to leave now?"/>
                    <br />
                    <br />
                    <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                        <font size="6"> ID : {window.batchid}</font>
                        &nbsp;&nbsp;&nbsp;&nbsp;
                        <Button variant="contained" color="default" onClick={this.getIdFile}>
                            Get ID
                        </Button>
                    </div>
                    <br />
                    <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                        <font size="4">Database : {window.databaseName}</font>
                    </div>
                    <br />
                    <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                        <font size="4">File name : {this.state.fileName}</font>
                    </div>
                    <br />
                    <br />
                    <br />
                    <br />
                    <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                        <CircularProgress size={175} />
                    </div>
                    <br />
                    <br />
                    <br />
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
                        <br />
                        <br />
                        <br />
                        <br />
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
                        </div>
    					<br />
    					<br />
                        <br />
                        <br />
                        <br />
                        <br />
                        <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
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
    	}
        
    }
}