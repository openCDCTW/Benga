import React from 'react';
import ReactDOM from 'react-dom';
import { withStyles } from '@material-ui/core/styles';
import SearchBar from 'material-ui-search-bar';
import QueryProfile from './QueryProfile.jsx';
import QueryDendrogram from './QueryDendrogram.jsx';

const styles = theme => ({
    displayCenter:{
        display:'flex',
        justifyContent:'center',
        alignItems:'center',
    },
})

class QuerybyID extends React.Component {

	constructor(props) {
		super(props);
		this.state = {};
	}

    queryProfile(){
        fetch('api/profiling/profile/' + this.state.queryProfile + '/', { method:'GET'})
        .then(function(response){
            if(response.status != 404){
                return response.json();
            }else{
                return response.status;
            }
        }).then(res => this.setState(state => ({ profile: res })));
    }

    queryDendrogram(){
        fetch('api/dendrogram/dendrogram/' + this.state.queryDendrogram + '/', { method:'GET'})
        .then(function(response){
            if(response.status != 404){
                return response.json();
            }else{
                return response.status;
            }
        }).then(res => this.setState(state => ({ dendrogram: res })));
    }

	render() {
		const { classes } = this.props;

    	return (
            <div style={{ marginTop: '70px'}}>
                <h2 className={classes.displayCenter}>Profile</h2>
                <SearchBar
                    onChange = {(value) => this.setState({ queryProfile: value })}
                    onRequestSearch={this.queryProfile.bind(this)}
                    placeholder = {"Input ID to get Profiles"}
                    style = {{
                        width: '70%',
                        margin: '0 auto',
                    }}
                />
                <QueryProfile profile={this.state.profile}/>
                <h2 className={classes.displayCenter}>Clustering Result</h2>
                <SearchBar
                    onChange = {(value) => this.setState({ queryDendrogram: value })}
                    onRequestSearch={this.queryDendrogram.bind(this)}
                    placeholder = {"Input ID to get clustering Result"}
                    style = {{
                        width: '70%',
                        margin: '0 auto',
                    }}
                />
                <QueryDendrogram dendrogram={this.state.dendrogram}/>
            </div>
    	)
    }
}

export default withStyles(styles)(QuerybyID)