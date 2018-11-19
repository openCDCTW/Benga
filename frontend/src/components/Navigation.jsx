import React from 'react';
import ReactDOM from 'react-dom';
import { Link } from 'react-router-dom';

export default class Navigation extends React.Component {

    render(){
        return (
        	<ul>
        	Upload  > 
        	<Link to="/profile_view" style={{textDecoration:'none'}}>
        		Profile  > 
        	</Link>
        	<Link to="/dendrogram_view">Dendrogram</Link>
        	</ul>
        );
    }
}