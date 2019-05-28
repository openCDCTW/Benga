import React from 'react';
import ReactDOM from 'react-dom';

export default class Header extends React.Component {


    render() {
        return (
            <header>
                <br />
                <div style={{ display:'flex', justifyContent:'center', 
                alignItems:'center'}}>
                    <a href="" target='_blank'>
                        <img src={require('./static/benga.png')} />
                    </a>
                </div>
                <a href="" style={{ textDecoration:'none' }} target='_blank'>
                	<font size="5" color="white">cgMLST@TAIWAN</font>
                </a>
                <br />
            </header>
        );
    }
}
