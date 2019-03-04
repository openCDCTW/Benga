import React from 'react';
import ReactDOM from 'react-dom';
import { Link,NavLink } from 'react-router-dom';
import Paper from '@material-ui/core/Paper';
import Tabs from '@material-ui/core/Tabs';
import Tab from '@material-ui/core/Tab';

export default class Navigation extends React.Component {

    constructor(props) {
        super(props);
        this.state = { value: 1 };
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
                    <Tab style={{ fontSize:'16px', textTransform:'none' }} disabled={window.tabSwitch} 
                    label="About" component={Link} to="/about" />
                    <Tab style={{ fontSize:'16px', textTransform:'none' }} disabled={window.tabSwitch} 
                    label="cgMLST Profiling" component={Link} to="/" />
                    <Tab style={{ fontSize:'16px',textTransform:'none' }} disabled={window.tabSwitch} 
                    label="Strain Tracking" component={Link} to="/tracking" />
                    <Tab style={{ fontSize:'16px',textTransform:'none' }} disabled={window.tabSwitch}
                    label="Search Database" component={Link} to="/tracking_search" />
                    <Tab style={{ fontSize:'16px',textTransform:'none' }} disabled={window.tabSwitch} 
                    label="Clustering" component={Link} to="/upload_profile" />
                    <Tab style={{ fontSize:'16px',textTransform:'none' }} disabled={window.tabSwitch} 
                    label="Example" component={Link} to="/demo" />
                    <Tab style={{ fontSize:'16px',textTransform:'none' }} disabled={window.tabSwitch}
                    label="Tutorial" component={Link} to="/tutorial" />
                </Tabs>
            </Paper>
        );
    }
}