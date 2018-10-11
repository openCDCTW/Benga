import React from 'react';
import ReactDOM from 'react-dom';

export default class Options_Scheme extends React.Component {


    render(){
        return (
        		<div>
        		<label>Scheme: &nbsp;&nbsp;</label>
        			<select name="choose_scheme" id="scheme">
        			<option value="95">95</option>
        			<option value="90">90</option>
        			<option value="85">85</option>
        			<option value="80">80</option>
        			<option value="50">50</option>
        			</select>
        		</div>
        );
    }
}