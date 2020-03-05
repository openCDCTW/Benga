import React from 'react';
import ReactDOM from 'react-dom';

export default class Footer extends React.Component {

	constructor(props) {
        super(props);
        window.year = new Date().getFullYear();
	};

    render() {
        return (
    		<footer style={{ marginTop:'40px', marginBottom: '40px',}}>
                <font style={{ fontSize: '14px'}}>
                    We recommend Google Chrome™ or Firefox® browser and set a screen resolution of minimum 1024 x 768 pixels for optimal performance.
                </font>
                <br />
                <font>Copyright (c) {window.year} Centers for Disease Control, MOHW, Taiwan (Taiwan CDC).</font>
      		</footer>
        );
    }
}