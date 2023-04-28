import { ComponentFixture, TestBed } from '@angular/core/testing';

import { ExpressionVisualizationComponent } from './expression-visualization.component';

describe('ExpressionVisualizationComponent', () => {
  let component: ExpressionVisualizationComponent;
  let fixture: ComponentFixture<ExpressionVisualizationComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      declarations: [ ExpressionVisualizationComponent ]
    })
    .compileComponents();

    fixture = TestBed.createComponent(ExpressionVisualizationComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
