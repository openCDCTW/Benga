import React from 'react';
import ReactDOM from 'react-dom';
import { withStyles } from '@material-ui/core/styles';
import { Link } from 'react-router-dom';
import { Prompt } from 'react-router';
import Button from '@material-ui/core/Button';
import ReplyIcon from '@material-ui/icons/Reply';
import DownloadIcon from '@material-ui/icons/CloudDownload';
import CircularProgress from '@material-ui/core/CircularProgress';
import download from 'downloadjs';

const styles = theme => ({
    displayCenter:{
        display:'flex',
        justifyContent:'center',
        alignItems:'center',
    },
})

class ProfilingResult extends React.Component {

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

        this.state = { fileName: fileName_str };
		this.queryProfile = this.queryProfile.bind(this);
	}

    queryProfile(){
		if(this.state.profile == undefined){
			fetch('api/profiling/profile/' + window.batchid, { method:'GET'})
			.then(response => response.json())
			.then(result => this.setState(state => ({
                profile: result.zip })));
		}else{
			clearInterval(this.interval);
		}
	}

    downloadJobID(){
        download(window.batchid,'BatchId.txt',"text/tab-separated-values");
    }

	componentDidMount(){
		this.interval = setInterval(this.queryProfile, 30000);
	}

	render() {
		const { classes } = this.props;
        if(this.state.profile == undefined){
    		return(
                <div style={{ marginTop: '100px', }}>
                    <Prompt
                        when={true}
                        message="You are leaving the page. Please save ID to get result. Are you sure to leave now?"/>
                    <div className={classes.displayCenter}>
                        <font size="6">Job ID : {window.batchid}</font>
                        <Button style={{ marginLeft: '20px', }} variant="contained" color="default" onClick={this.downloadJobID}>
                            Get job ID
                        </Button>
                    </div>
                    <br />
                    <div className={classes.displayCenter}>
                        <font size="4">Database : {window.database}</font>
                    </div>
                    <br />
                    <div className={classes.displayCenter}>
                        <font size="4">File name : {this.state.fileName}</font>
                    </div>
                    <div style={{ marginTop: '100px', marginBottom: '100px' }} className={classes.displayCenter}>
                        <CircularProgress size={175} />
                    </div>
    			</div>
        	);
        }else{
    		return (
                <div>
                    <Prompt
                        when={true}
                        message="You are leaving the page. Please save results, or it will lose. Are you sure to leave now?"/>
                    <div style={{ marginTop: '100px' }} className={classes.displayCenter}>
                        <font size="6">Job ID : {window.batchid}</font>
                        <Button style={{ marginLeft: '20px' }} variant="contained" color="default" onClick={this.downloadJobID}>
                            Get job ID
                        </Button>
                    </div>
                    <div className={classes.displayCenter}>
                        <a download href={this.state.profile}
                         style={{ textDecoration:'none', margin: '70px', }}>
                            <Button variant="contained" color="default">
                            Download profiles (.zip)
                            &nbsp;&nbsp;
                            <DownloadIcon />
                            </Button>
                        </a>
                    </div>
                    <div className={classes.displayCenter}>
                        <Link to="/cgMLST" style={{ textDecoration:'none' }}>
                            <Button variant="contained" color="default">
                                <ReplyIcon />
                                &nbsp;&nbsp;
                                Back
                            </Button>
                        </Link>
                    </div>
                </div>
    		);
        }

    }
}

export default withStyles(styles)(ProfilingResult)