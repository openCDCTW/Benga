import React from 'react';
import ReactDOM from 'react-dom';
import { BrowserRouter as Router, Route } from "react-router-dom";
import Header from './Header.jsx';
import Footer from './Footer.jsx';
import Upload from './Upload.jsx';
import Output from './Output.jsx'

// Render

class Main extends React.Component {


    render() {
        return(
            // Router：在最外層的路由組件
            <Router>
                {/* Router底下只能有一個child */}
                <div>
                    {/* Route：路由的配對組件 */}
                    {/* path是要與此路由配對的路徑、exact是指是否要路徑完成匹配、component是指render的組件 */}

                    <Header />
                    <Route path="/" exact component={Upload} />
                    <Route path="/profiling" component={Output} />
                    <Footer />
                </div>
            </Router>
        );
    }
}

ReactDOM.render(<Main /> ,document.getElementById('root'));
