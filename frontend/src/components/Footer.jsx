import React from 'react';
import ReactDOM from 'react-dom';

export default class Footer extends React.Component {

	constructor(props) {
            super(props);
	};

    render() {
        return (
    		<footer>
                <br />
                <br />
                <br />
        		<div>
                    <br />
                    <font style={{ fontSize: '14px'}}>
                        We recommend Google Chrome™ or Firefox® browser and set a screen resolution of minimum 1024 x 768 pixels for optimal performance.
                    </font>
                    <br />
                    <br />
                    <font>Copyright (c) 2020 Centers for Disease Control, MOHW, Taiwan (Taiwan CDC).</font>
                </div>
                <br />
                <br />
                <br />
      		</footer>
        );
    }
}
