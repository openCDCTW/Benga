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

export default class Profile_view extends React.Component {

	constructor(props) {
		super(props);
		this.state = {};
		this.query_profile_result = this.query_profile_result.bind(this);
	};

	query_profile_result(){

		if(this.state.profile_result_Url == undefined){
			fetch('api/profiling/profile/' + window.batchid, { method:'GET'})
			.then(response => response.json())
			.then(result => this.setState(state => ({profile_result_Url: result.file})));
		}else{
			clearInterval(this.interval);
		}

	}

	componentDidMount(){
		this.query_profile_result();
		this.interval = setInterval(this.query_profile_result, 60000);
	}

    render() {

    	if(this.state.profile_result_Url == undefined){

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
                    <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                        <font> Please hold on ... </font>
                    </div>
                    <br />
                    <div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
                        <img src="https://svgshare.com/i/9N5.svg" />
                    </div>
                    <br />
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
                            <a download href={this.state.profile_result_Url} 
                             style={{ textDecoration:'none' }}>
                                <Button variant="contained" color="default">
                                Download all profiles 
                                &nbsp;&nbsp;
                                <DownloadIcon />
                                </Button>
                            </a>
                        </div>
    					<br />
    					<br />
                        <Link to="/" style={{ textDecoration:'none' }}>
                            <Button variant="contained" color="default">
                                <ReplyIcon />
                                &nbsp;&nbsp;
                                Back
                            </Button>
                        </Link>
    				</div>
        	);
    	}
        
    }
}
