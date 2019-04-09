import React from 'react';
import ReactDOM from 'react-dom';

export default class Footer extends React.Component {


    render() {
        return (
    		<footer>
        		<div>
                    <font>Copyright (c) 2019 Centers for Disease Control, MOHW, Taiwan (Taiwan CDC).</font>
                    <br />
                    <br />
                    <font style={{ fontSize: '14px'}}>
                    	建議您使用Google Chrome、Firefox瀏覽器，並搭配 1024 x 768 以上之螢幕解析度，以獲得最佳瀏覽體驗。
                    </font>
                    <br />
                    <font style={{ fontSize: '14px'}}>
                        We recommend using Google Chrome, Firefox and with a resolution of 1024 x 768 (or above) for the best browsing experience.
                    </font>
                </div>
      		</footer>
        );
    }
}