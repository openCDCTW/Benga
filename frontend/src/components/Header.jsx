import React from 'react';
import ReactDOM from 'react-dom';

export default class Header extends React.Component {


    render() {
        return (
            <header>
                <br />
                <a href="http://localhost:8000" style={{ textDecoration:'none' }} target='_blank'>
                	<font size="5" color="white">BENGA cgMLST@TAIWAN</font>
                </a>
                <br />
            </header>
        );
    }
}
