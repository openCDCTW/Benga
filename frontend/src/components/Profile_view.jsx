import React from 'react';
import ReactDOM from 'react-dom';
import { Link } from 'react-router-dom';

export default class Profile_view extends React.Component {

	constructor(props) {
		super(props);
		this.state = {};
		this.query_profile_result = this.query_profile_result.bind(this);
	};

	query_profile_result(){

		if(this.state.profile_result_Url == undefined){
			fetch('api/profiling/profile/'+window.batchid, {method:'GET'})
			.then(response => response.json())
			.then(result => this.setState(state=>({profile_result_Url:result.file})));
		}else{
			clearInterval(this.interval);
		}

	}

	componentDidMount(){
		this.query_profile_result();
		this.interval = setInterval(this.query_profile_result, 30000);
	}

    render() {

    	if(this.state.profile_result_Url == undefined){

    		return(<img src="http://i.imgur.com/1J9vdDD.gif" alt="Please waiting..." />);
    	
    	}else{
    		return (
    				<div id="url">
    				<br />
    				<a download href={this.state.profile_result_Url}>Download all profiles</a>
    				<br />
    				</div>
        	);
    	}
        
    }
}
