import React from 'react';
import ReactDOM from 'react-dom';

export default class Header extends React.Component {


    render() {
        return (
            <header style={{ backgroundColor: '#2828ff' }}>
                <br />
                <div style={{ display:'flex', justifyContent:'center', 
                alignItems:'center'}}>
                    <a href="http://rdvd.cdc.gov.tw/cgMLST/non-release" target='_blank'>
                        <img src={require('../components/static/benga.png')} />
                    </a>
                </div>
                <a href="http://rdvd.cdc.gov.tw/cgMLST/non-release" style={{ textDecoration:'none' }} target='_blank'>
                        <font size="5" color="white">cgMLST@TAIWAN</font>
                </a>
                <br />
            </header>
        );
    }
}
