import React from 'react';
import ReactDOM from 'react-dom';
import PropTypes from 'prop-types';
import { withStyles } from '@material-ui/core/styles';
import InputLabel from '@material-ui/core/InputLabel';
import MenuItem from '@material-ui/core/MenuItem';
import FormHelperText from '@material-ui/core/FormHelperText';
import FormControl from '@material-ui/core/FormControl';
import Select from '@material-ui/core/Select';

const styles = theme => ({
  root: {
    display: 'flex',
    flexWrap: 'wrap',
  },
  formControl: {
    marginLeft: theme.spacing.unit,
    minWidth: 240,
  },
  selectEmpty: {
    marginTop: theme.spacing.unit * 2,
  },
});

class Options extends React.Component {

    constructor(props) {
        super(props);

        this.state = {
            database: 'Vibrio_cholerae 95',
            labelWidth: 0,
        };

        this.handleChange = this.handleChange.bind(this);
    }
    
    handleChange(event){
	let db_info = event.target.value.split(" ")
        this.setState(state => ({ [event.target.name]: event.target.value}));
	
	window.databaseName = db_info[0];
	window.occurrence = db_info[1];
    };

    render(){

        const { classes } = this.props;
        
        return (
          <form className={classes.root} autoComplete="off">
            <FormControl required className={classes.formControl} disabled={this.props.switch}>
              <InputLabel htmlFor="database-required">Species</InputLabel>
                <Select
                  value={this.state.database}
                  onChange={this.handleChange}
                  name="database"
                  inputProps={{
                    id: 'database-required',
                  }}
                  className={classes.selectEmpty}
                  >
                  <MenuItem value={'Campylobacter_jejunii/coli'} disabled>Campylobacter coli/jejunii</MenuItem>
                  <MenuItem value={'Cronobacter_sakazakii 95'}>Cronobacter sakazakii</MenuItem>
                  <MenuItem value={'Escherichia_coli/Shigella'} disabled>Escherichia coli/Shigella</MenuItem>
                  <MenuItem value={'Listeria_monocytogenes 95'}>Listeria monocytogenes</MenuItem>
                  <MenuItem value={'Salmonella_enterica 94'}>Salmonella enterica</MenuItem>
                  <MenuItem value={'Vibrio_cholerae 95'}>Vibrio cholerae</MenuItem>
                  <MenuItem value={'Vibrio_parahaemolyticus'} disabled>Vibrio parahaemolyticus</MenuItem>
                </Select>
              <FormHelperText>Required</FormHelperText>
            </FormControl>
          </form>
        );
    }
}

Options.propTypes = {
  classes: PropTypes.object.isRequired,
};

export default withStyles(styles)(Options);

// style={{ color:'#d0d0d0'}}
