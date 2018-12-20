import React from 'react';
import ReactDOM from 'react-dom';
import { Link,NavLink } from 'react-router-dom';
import Paper from '@material-ui/core/Paper';
import Tabs from '@material-ui/core/Tabs';
import Tab from '@material-ui/core/Tab';

export default class Navigation extends React.Component {

    constructor(props) {
        super(props);
    };

    render(){
        return (
        	<Paper square>
                <Tabs value={this.props.value} indicatorColor="primary" textColor="primary" centered>
                    <Tab label="Profile" component={Link} to="/" />
                    <Tab label="Dendrogram" component={Link} to="/upload_profile" />
                    <Tab label="Example" component={Link} to="/demo" />
                </Tabs>
            </Paper>
        );
    }
}

