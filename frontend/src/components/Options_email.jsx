import React from 'react';
import ReactDOM from 'react-dom';

export default class Email extends React.Component {


    render(){
        return (
        	<div>
            <label>Email:&nbsp;&nbsp;&nbsp; </label>
            <input id="email" name="email" type="text" placeholder="Input your Email here"/>
            </div>
        );
    }

    
    
}
