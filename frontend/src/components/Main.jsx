import React from 'react';
import ReactDOM from 'react-dom';
import { BrowserRouter as Router, Route } from "react-router-dom";
import Header from './Header.jsx';
import Navigation from './Navigation.jsx';
import Footer from './Footer.jsx';
import Upload from './Upload.jsx';
import Profile_view from './Profile_view.jsx';
import Dendrogram_view from './Dendrogram_view.jsx'

class Main extends React.Component {

    render() {
        return(
            <Router>
                <div>
                    <Header />
                    <Navigation />
                    <Route path="/" exact component={Upload} />
                    <Route path="/profile_view" component={Profile_view} />
                    <Route path="/dendrogram_view" component={Dendrogram_view} />
                    <Footer />
                </div>
            </Router>
        );
    }
}

ReactDOM.render(<Main /> ,document.getElementById('root'));
