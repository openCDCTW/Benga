import React from 'react';
import ReactDOM from 'react-dom';
import { Link } from 'react-router-dom';
import Paper from '@material-ui/core/Paper';
import Tabs from '@material-ui/core/Tabs';
import Tab from '@material-ui/core/Tab';
import Button from '@material-ui/core/Button';
import { withStyle } from '@material-ui/core/styles';
import ReplyIcon from '@material-ui/icons/Reply';

export default class Tracking_result extends React.Component {

	constructor(props) {
		super(props);
		this.state = {};
		this.query_track_result = this.query_track_result.bind(this);
	};

	query_track_result(){


	}

	componentDidMount(){
		this.query_track_result();
		this.interval = setInterval(this.query_track_result, 10000);
	}



    render() {

    	if(this.state.tracking_result == undefined){

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
						<Tabs value={false} centered>
							<Tab label=" " disabled/>
						</Tabs>
					</Paper>
					<br />
					<div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
					</div>
					<br />
					<br />
					<br />
					<br />
					<br />
					<br />
					<div style={{ display:'flex', justifyContent:'center', alignItems:'center'}}>
						<Link to="/tracking" style={{ textDecoration:'none' }}>
							<Button variant="contained" color="default">
								<ReplyIcon />
								&nbsp;&nbsp;
								Back
							</Button>
						</Link>
					</div>
					<br />
				</div>
			);
		}
    }
}
