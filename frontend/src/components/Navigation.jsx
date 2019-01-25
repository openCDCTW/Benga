import React from 'react';
import ReactDOM from 'react-dom';
import { Link,NavLink } from 'react-router-dom';
import Paper from '@material-ui/core/Paper';
import Tabs from '@material-ui/core/Tabs';
import Tab from '@material-ui/core/Tab';

export default class Navigation extends React.Component {

    constructor(props) {
        super(props);
        this.state = { value: 0 };
        window.tabSwitch = false;
    };

    handleChange(event, value){
        this.setState({ value });
    };

    render(){
        const { classes } = this.props;
        const { value } = this.state.value;

        return (
        	<Paper square>
                <Tabs value={this.state.value} indicatorColor="primary" 
                textColor="primary" scrollButtons="auto" scrollable={true} 
                onChange={this.handleChange.bind(this)} disabled={window.tabSwitch}>
                    <Tab disabled={window.tabSwitch} label="Profile" component={Link} to="/" />
                    <Tab disabled={window.tabSwitch} label="Dendrogram" component={Link} to="/upload_profile" />
                    <Tab disabled={window.tabSwitch} label="Tracking" component={Link} to="/tracking" />
                    <Tab disabled={window.tabSwitch} label="Example" component={Link} to="/demo" />
                    <Tab disabled={window.tabSwitch} label="Tutorial" component={Link} to="/tutorial" />
                </Tabs>
            </Paper>
        );
    }
}

// <Tab disabled={window.tabSwitch} label="Search" component={Link} to="/tracking_search" />